library(tidyverse)

### Read cover data 1KB windows
#####

depth <- read_csv("../data/mean_read_counts_1kb.csv")%>%
  separate(id, c("scaffold", "start"),sep="_", remove=F)%>%
  mutate(start = as.numeric(start))%>%
  rename(raw_f = 'f', raw_m = 'm')%>%
  mutate(f = 1e6 * raw_f / sum(raw_f), m = 1e6 * raw_m/ sum(raw_m) )

depth_10k <- depth %>% 
  group_by(scaffold)%>%
  mutate(win = floor(start/10000)+1)%>%
  group_by(scaffold, win)%>%
  summarize(f = mean(f, na.rm=T), m=mean(m, na.rm=T), start = min(start))%>%
  ungroup()%>%
  mutate(win=NULL)
  
win_1k <- depth %>% transmute(id, scaffold, start)
win_10k <- win_1k %>% filter(start%%10000 == 1)
win_20k <- win_1k %>% filter(start%%20000 == 1)

#####
### PopGen Stats 10KB windows
#####

# Fst and Dxy btw males and females in two pops
raw_MF670 <- read_csv(gzfile("../data/FstDxy_cal/670_sex-diff.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  rename(pi_M670 = "pi_popB", pi_F670 =  "pi_popA", dxy_MF670 = "dxy_popA_popB", fst_MF670 = "Fst_popA_popB" )

raw_MF716 <- read_csv(gzfile("../data/FstDxy_cal/716_sex-diff.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  rename(pi_M716 = "pi_popB", pi_F716 =  "pi_popA", dxy_MF716 = "dxy_popA_popB", fst_MF716 = "Fst_popA_popB" )

# Fst and Dxy btw two pops
raw_popdiff <- read_csv(gzfile("../data/FstDxy_cal/670_vs_716.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  rename(pi_670 = "pi_popA", pi_716 =  "pi_popB", popDxy = "dxy_popA_popB", popFst = "Fst_popA_popB" )

# merge the three dataframes
pgstats <- left_join(raw_MF670, raw_MF716) %>% left_join(raw_popdiff) %>%
  mutate(popFst = if_else(popFst < 0, 0, popFst),
         fst_MF670 = if_else(fst_MF670 < 0, 0, fst_MF670),
         fst_MF716 = if_else(fst_MF716 < 0, 0, fst_MF716))%>%
  mutate(sexDxy = (dxy_MF670 + dxy_MF716)/2 )%>%
  mutate(sexFst = (fst_MF670 + fst_MF716)/2)%>%
  mutate( piF = 0.5*(pi_F670 + pi_F716),
          piM = 0.5*(pi_M670 + pi_M716),
          pipop =0.5*(pi_670 + pi_716) )

pgstats_win  <-  win_10k %>% left_join(pgstats, by=c('id', 'scaffold', 'start')) %>%
  mutate(end=start+1e4)%>%
  replace_na(list(sexFst=0, piF=0, piM=0, sexDxy=0))

#####
### SSR count data 10KB windows
#####
fsr_raw <- read_delim("../data/SSR_mapped_data/filtered_FSR_10kb.txt",delim="\t")%>%
  rename(fsr = 'read_mapped', fsr_breadth = 'breadth')%>%
  mutate(covered_sites = NULL, total_sites = NULL)
msr_raw <- read_delim("../data/SSR_mapped_data/filtered_MSR_10kb.txt",delim="\t")%>%
  rename(msr = 'read_mapped', msr_breadth = 'breadth')%>%
  mutate(covered_sites = NULL, total_sites = NULL)

ssr_clean <- full_join(fsr_raw, msr_raw, by=c('scaffold', 'start', 'end'))%>%
  gather(sex, value = read_count, one_of('fsr', 'msr'))%>%
  mutate(breadth=map2_dbl(msr_breadth, fsr_breadth, max, na.rm=T), msr_breadth=NULL, fsr_breadth=NULL)%>%
  na.omit()%>%
  mutate(start= start+1)%>%
  unite("id", scaffold, start,sep="_", remove=F)

ssr_10kb <- win_10k %>% 
  left_join(ssr_clean%>%transmute(id,sex,read_count,breadth), by='id')%>%
  replace_na(list(read_count=0, breadth=0))

#####
### Gene annotation
#####
raw_annot <- read_delim("../data/Genome_annotation/ahrd_output.txt", skip = 2,  delim="\t")%>%
  separate(`Protein-Accession`,"genename", sep="-", extra = "drop", remove=F)%>%
  transmute(protein =`Protein-Accession`, genename, bhit = `Blast-Hit-Accession`, description = `Human-Readable-Description` )

raw_locs <- read_delim("../data/Genome_annotation/sapp.aed.0.5.gff", 
                       delim="\t", 
                       col_names = c("scaffold", "x1",'x2','start', 'end', 'q','strand', 'x3', 'annot'))%>%
  filter(x2 =='gene')%>%
  separate(annot, c("id"), sep=";",extra = "drop")%>%
  separate(id, c("idtag", "genename"), sep="=",extra = "drop")%>%
  transmute(genename, scaffold, start, end)

annot <- raw_locs%>% left_join(raw_annot, by='genename')%>%
  select(one_of("scaffold", "genename", "bhit", 'description', 'start', 'end'))%>%
  distinct()

#####
### Expression data
#####
raw_flor <- read_csv("../data/Expr_data/allGene_Expr_flower.csv")%>%
  rename(genename = 'Geneid')

raw_buds <- read_csv("../data/Expr_data/allGene_Expr_bud.csv")%>%
  rename(genename = 'Geneid')

degs_f <- raw_flor%>%
  transmute(genename, logFC, cpm=logCPM, pval=PValue)%>%
  left_join(raw_locs, by='genename')%>%mutate(tissue='f')

degs_b <- raw_buds%>%
  transmute(genename, logFC, cpm=logCPM, pval=PValue)%>%
  left_join(raw_locs, by='genename')%>%mutate(tissue='b')

express_threshold <- 0.1

degs_f_filtered <- degs_f %>% filter(cpm > quantile(cpm, c(express_threshold)))
degs_b_filtered <- degs_b %>% filter(cpm > quantile(cpm, c(express_threshold)))

degs <- bind_rows(degs_f_filtered, degs_b_filtered)

degs_10k <- degs %>% 
  mutate(win_start = floor(start/10000)*10000+1)%>%
  unite("id", scaffold, win_start, sep="_", remove=F)%>%
  group_by(id)%>%
  summarise(avg_exp = mean(cpm), genecount = n(), sigcount = sum(pval < 1e-2))%>%
  right_join(win_10k, by='id')%>%
  replace_na(list(avg_exp=0, genecount=0, sigcount=0))


#####
### Lyco orthologs
#####

itag <- read_delim("../data/ITAG3.2_brief_mappings.txt", delim="\t", col_names=c("start", "end", "gene"))

blastraw <- tibble( filename =list.files("blastOut/"), fullname = str_c("blastOut/",filename, sep="")) %>%
  mutate(content = map(fullname, read_delim, delim="\t", col_names = FALSE))%>%
  unnest()%>%
  separate(col = filename, into = c("scaffold", "ext"),sep = "_")%>%
  transmute(scaffold, sappname = X1, solname = X2)

orthologs <- blastraw %>%
  separate(solname, "lycochr", sep="g", remove = F, extra='drop')%>%
  separate(solname, "lyconame", sep="\\.", extra='drop')%>%
  mutate(lycochr = str_replace(lycochr, 'Solyc',""))%>%
  separate(sappname, "genename", sep="-" , extra='drop' )%>%
  distinct()%>% 
  mutate(start = annot$start[match(genename, annot$genename)])%>%
  mutate(lycopos = itag$start[match(lyconame, itag$gene)])

lycochr_infer <- orthologs %>%
  group_by(scaffold, lycochr)%>%
  summarize(ngenes = n())%>%
  mutate(frac = ngenes/sum(ngenes))%>%
  arrange(scaffold, desc(ngenes))%>%
  filter(row_number()<2)%>%
  ungroup()

syntenymap <- orthologs %>%
  mutate(chrom_inferred = lycochr_infer$lycochr[match(scaffold, lycochr_infer$scaffold)],
         pos_inferred = if_else(chrom_inferred == lycochr, lycopos, 0))

total_win_df <- full_join(pgstats_win, ssr_10kb, by=c("id", "scaffold", "start")) 


#####
## Functions
#####

genegrab_fn <- function(scaf, ws, we, df, windowed = F){
  if(windowed){
    df%>%
      filter(scaffold %in% scaf)%>%
      filter(start >= ws & start < we)
  }
  else{
    df%>%
      filter(scaffold %in% scaf)%>%
      filter( (start >= ws & start < we) | (end >= ws & end < we))
  }
}
ggrab_df  <- function(df, regions, windowed=F){
  regions %>%
    mutate(pgs = pmap(list(scaf=scaffold, ws=winstart, we=winend), genegrab_fn, df = df, windowed=windowed))%>%
    ungroup()%>%
    unnest()
}

myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

#####
### Filtering SSRs
#####
ssr <- ssr_clean %>%
  mutate(sexFst = pgstats$sexFst[match(id,pgstats$id)])%>%
  mutate(sexDxy = pgstats$sexDxy[match(id,pgstats$id)])%>%
  mutate(popDxy = pgstats$popDxy[match(id,pgstats$id)])%>%
  mutate(popFst = pgstats$popFst[match(id,pgstats$id)])%>%
  na.omit()

ssr_scaffs <- ssr%>%group_by(scaffold)%>%summarise(n=n())

#How many SSRs? 153
length(ssr$id)
#How many scaffolds carry SSRs? 65
length(ssr_scaffs$scaffold)

ssr_top <- ssr %>%
  filter(breadth > 0.1)
ssr_top_scaffs<- ssr_top%>%
  group_by(scaffold,sex)%>%summarise(n=n())%>%arrange(desc(n))

#How many top SSRs? 18
length(ssr_top$id)
#How many scaffolds carry top SSRs? 10
length(ssr_top_scaffs$scaffold)
#write_csv(arrange(ssr_top_scaffs, desc(n) )%>%select(scaffold), "top_SSR_scaffolds.csv")

div_threshold <- 0.9
quant_sexDxy <- quantile(pgstats$sexDxy,c(div_threshold), na.rm=T)
quant_sexFst<- quantile(pgstats$sexFst,c(div_threshold), na.rm=T)

ssr_div <- ssr_top %>%
  filter(sexFst > popFst) %>%  
  filter(sexDxy > popDxy) %>%  
  filter(sexDxy > quant_sexDxy)%>%
  filter(sexFst > quant_sexFst)%>%
  arrange(scaffold)

ssr_div_scaffs <- ssr_div %>%
  group_by(scaffold,sex)%>%summarise(n=n())%>%arrange(desc(n))

#How many top SSRs with Fst? 
length(ssr_div$id)
#How many scaffolds carry top SSRs with Fst? 
length(ssr_div_scaffs$scaffold)



#####
## export around windows
#####

ssk_export <- ssr_div %>% 
  mutate(scaffold = paste("scf", str_remove(scaffold, "scf71800000"), sep="_" ))%>%
  transmute(sex, winname = paste(scaffold, ":", start,"-", end, sep=""), breadth, sexFst, sexDxy)

#write_csv(ssk_export, "top_ssk_stats.csv")

ssk_regions <- ssr_div %>%
  transmute(scaffold, winstart = start - 2e4, winend =  end + 2e4)

inside_genes <- ggrab_df(annot%>%transmute(scaffold, start, end, genename), ssk_regions)%>%
  transmute(scaffold, start, genename, type="in")%>%
  distinct()%>%
  left_join(degs%>%select(genename, cpm, tissue), by=c('genename'))%>%
  spread(tissue, cpm)%>%
  mutate(`<NA>`=NULL)%>%
  mutate(scaffold = paste("scf", str_remove(scaffold, "scf71800000"), sep="_" ))%>%
  mutate(description = annot$description[match(genename, annot$genename)],
         blast_hit = orthologs$lyconame[match(genename, orthologs$genename)])%>%
  arrange(scaffold, start)

#write_csv(inside_genes, "genes_20KB_around_SSR.csv")


ssk_wider <- ssr_div %>%
  transmute(scaffold, winstart = start - 1e5, winend =  end + 1e5)

around_genes <- ggrab_df(degs, ssk_wider) %>% 
  filter(pval <0.01 | cpm > quantile(degs$cpm, c(0.25))) %>%
  myspread(tissue, c(logFC, cpm, pval))%>%
  transmute(scaffold, start, genename, b_cpm, b_logFC, b_pval, f_cpm, f_logFC, f_pval)%>%
  mutate(scaffold = paste("scf", str_remove(scaffold, "scf71800000"), sep="_" ))%>%
  distinct()%>%
  mutate(description = annot$description[match(genename, annot$genename)],
         blast_hit = orthologs$lyconame[match(genename, orthologs$genename)])%>%
  arrange(scaffold, start)

#write_csv(around_genes, "genes_100KB_around_SSR.csv")

### mu value for tomato: 1.3 x10-8 (From Rice; Ma & Bennetzen 2004)
