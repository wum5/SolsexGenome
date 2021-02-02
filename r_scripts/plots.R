library(patchwork)
#####
###Plots
#####

window_theme_mine <- function(p, regions, placement=0){
  
  p$layers <- c(geom_rect(data=regions, aes(xmin=start, xmax=end, ymin=Inf, ymax=-Inf), fill='grey85'), p$layers)
  out <- p+
    facet_grid(.~cleanscaff, scales = "free_x", space='free')+
    theme_classic()+
    theme(legend.position = 'none',                        
          strip.background = element_blank(), 
          panel.border = element_rect(fill=NA))
  
  
  notop <- theme(strip.text = element_blank())
  nobottom <- theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) #
  if(placement==1){
    out + nobottom
  }
  else if(placement==2){
    out+ 
      notop+
      scale_x_continuous(labels = function(x)x/1000)+
      xlab("Position (Kb)")
  }
  else{out + nobottom + notop}  
}

#mean(total_win_df$pipop, na.rm=T)
plotall_ssk <- function(ssr_in, pg_df = total_win_df, expr_df = degs, annot_df = annot, depth_df=depth, window_width = 5e4+1, pimean = 0.002800299, pb_df =pbssx){
  ssr_local <- ssr_in %>%
    mutate(sex = factor(sex, levels = c('fsr', 'msr'), labels=c("Female region", "Male region")))%>%
    mutate(cleanscaff = paste("scf", str_remove(scaffold, "scf71800000"), sep="_" ))
  
  closeup_region <- ssr_local %>%
    group_by(sex, scaffold, cleanscaff)%>%
    summarize(winstart = max(0, min(start) - window_width),
              winend =  max(start) + 1e4 + window_width)
  
  cut_wins <- ggrab_df(pg_df, closeup_region, T) 
  cut_depth <- ggrab_df(depth_df, closeup_region, T)
  cut_annot <- ggrab_df(annot_df, closeup_region, T) 
  cut_degs <- ggrab_df(expr_df, closeup_region, T) 
  cut_pbs <- ggrab_df(pb_df, closeup_region, T) 
  sigdegs <- sum(cut_degs$pval < 0.01)
  
  # piplot <- ggplot(cut_wins)+
  #   geom_step(aes(x=start,y = piM), color='blue')+
  #   geom_step(aes(x=start,y = piF), color='red')+
  #   geom_hline(yintercept = quantile(pg_df$pipop, probs = c(0.5), na.rm=T), lty=2)+
  #   geom_hline(yintercept = quantile(pg_df$pipop, probs = c(0.9), na.rm=T), lty=3)+
  #   labs(y="Heterozygosity")+
  #   scale_y_continuous(breaks=c(0.005,0.015))+
  #   expand_limits(y=0)

  piplot <-ggplot(cut_wins)+
    geom_hline(yintercept = 0)+
    geom_step(aes(x=start,y =0), color=NA)+
    geom_point(aes(x=start+5e3, y = (piM-piF)/pimean, color= (piM-piF)/pimean < 0), alpha=0.5)+
    geom_segment(aes(x=start+5e3, xend = start+5e3, y = 0, yend = (piM-piF)/pimean, color= (piM-piF)/pimean < 0), alpha=0.5)+
    labs(y=expression("M-F "*Delta*pi))+
    scale_color_manual(values=c("darkblue", "red"))
    expand_limits(y=0)
  
  fstplot<- ggplot(cut_wins)+
    geom_point(aes(x=start+5e3, y=0), fill=NA, color=NA)+
    geom_step(aes(x=start,y =sexFst), color='purple')+
    geom_hline(yintercept = quantile(pg_df$sexFst, probs = c(0.5), na.rm=T), lty=2)+
    geom_hline(yintercept = quantile(pg_df$sexFst, probs = c(0.9), na.rm=T), lty=3)+
    labs(y=expression("M-F "*italic(F)["ST"]))+
    scale_y_continuous(breaks=c(0, 0.2, 0.4))+
    expand_limits(y=0)
 
  dxyplot <- ggplot(cut_wins)+
    geom_point(aes(x=start+5e3, y=0), fill=NA, color=NA)+
    geom_step(aes(x=start,y = sexDxy), color='darkgreen', alpha=0.5)+
    #geom_step(aes(x=start,y = dxy_MF670), color='blue', alpha=0.5)+
    #geom_step(aes(x=start,y = popDxy), color='black', alpha=0.5)+
    geom_hline(yintercept = quantile(pg_df$sexDxy, probs = c(0.5), na.rm=T), lty=2)+
    geom_hline(yintercept = quantile(pg_df$sexDxy, probs = c(0.9), na.rm=T), lty=3)+
    labs(y=expression("M-F "*italic(d)["XY"]))+
    scale_y_continuous(breaks=c(0, 0.005,0.015))+
    expand_limits(y=0)
  
  depthplot<- ggplot(cut_depth)+  
    geom_point(data=cut_wins, aes(x=start+5e3, y=0), fill=NA, color=NA)+
    geom_step(aes(x=start, y =f), alpha=0.5, color='red')+
    geom_step(aes(x=start, y =m), alpha=0.5, color='blue')+
    ylab("Coverage")+
    coord_cartesian(ylim=c(0, 2.3))+
    scale_y_continuous(breaks = c(0,1,2))
  
  degplot <- ggplot(cut_degs)+
    geom_point(data=cut_wins, aes(x=start+5e3, y=0), fill=NA, color=NA)+
    geom_segment(data=cut_pbs, aes(x=winstart, xend=winend, y=0,yend=0))+
    geom_step(data=cut_wins, aes(x=start, y=breadth))+
    geom_point(data=cut_annot, aes(x=start, y=0.75, fill = genename %in% cut_degs$genename), shape= 22, alpha=0.5)+
    scale_fill_manual(values = c(NA,'black'))+
    ylab("SSR breadth")
  if(sigdegs >0){ degplot<- degplot + geom_point(data=filter(cut_degs,pval <0.01), aes(x=start, y=0.75), color='red', fill='red', shape=22, alpha=0.5)  }
  p1 <- window_theme_mine(degplot, ssr_local, 1)
  p2 <- window_theme_mine(depthplot, ssr_local, 0)
  p3 <- window_theme_mine(piplot, ssr_local, 0)
  p4 <- window_theme_mine(fstplot, ssr_local, 0)
  p5 <- window_theme_mine(dxyplot, ssr_local, 2)
  p1+p2+p3+p5+ patchwork::plot_layout(ncol=1, heights = c(0.5,0.25, 0.5,0.5))
}

ssr_tmp <- ssr_clean %>% arrange(desc(breadth))%>% filter(row_number()<11)

fivepanelplot <- plotall_ssk(ssr_tmp, window_width = 1e5+1)

ggsave("hd5plot.pdf",fivepanelplot, width=8, height=5, useDingbats=FALSE)

#plotall_ssk(ssr_div,  window_width = 5e6+1)
#####

plot_chr12 <- function(ssr_in, pg_df = total_win_df, expr_df = degs, annot_df = annot, depth_df=depth, synt_df = syntenymap, window_width = 5e4+1, pimean = 0.002800299){
  ssr_local <- ssr_in %>%
    mutate(sex = factor(sex, levels = c('fsr', 'msr'), labels=c("Female region", "Male region")))%>%
    mutate(cleanscaff = paste("scf", str_remove(scaffold, "scf71800000"), sep="_" ))
  
  closeup_region <- ssr_local %>%
    group_by(sex, scaffold, cleanscaff)%>%
    summarize(winstart = max(0, min(start) - window_width),
              winend =  max(start) + 1e4 + window_width)
  
  cut_wins <- ggrab_df(pg_df, closeup_region, T) 
  cut_depth <- ggrab_df(depth_df, closeup_region, T)
  cut_annot <- ggrab_df(annot_df, closeup_region, T) 
  cut_degs <- ggrab_df(expr_df, closeup_region, T) 
  cut_ortho <- ggrab_df(synt_df, closeup_region, T) 

  oplot <-ggplot(cut_ortho)+
    geom_point(aes(x=start, y=pos_inferred/1e6, shape= pos_inferred < 2e6))+
    scale_y_continuous(limits=c(2, 4.5), oob = scales::squish)+
    scale_shape_manual(values = c(16, 4))+
    labs(y=expression(paste("Position in ",italic("S. lycopersicum"),"(Mb)")))
  # 
  # fstplot<- ggplot(cut_wins)+
  #   geom_step(aes(x=start,y =sexFst), color='purple')+
  #   geom_hline(yintercept = quantile(pg_df$sexFst, probs = c(0.5), na.rm=T), lty=2)+
  #   geom_hline(yintercept = quantile(pg_df$sexFst, probs = c(0.9), na.rm=T), lty=3)+
  #   labs(y=expression("M-F "*italic(F)["ST"]))+
  #   scale_y_continuous(breaks=c(0.2, 0.4))+
  #   expand_limits(y=0)
  
  degplot <- ggplot(cut_wins)+
    geom_step(aes(x=start, y=breadth))+
    scale_shape_manual(values = c(0,15))+
    ylab("SSR breadth")
  
  p1 <- window_theme_mine(oplot, ssr_local, 1)
  p3 <- window_theme_mine(degplot, ssr_local, 2)
 # p4 <- window_theme_mine(fstplot, ssr_local, 2)
  p1+p3 + patchwork::plot_layout(ncol=1, heights = c(1,0.5))
}

scaffs_12 <- syntenymap %>% filter(chrom_inferred ==12)%>%pull(scaffold)%>%unique()

ssr_12 <- ssr_div %>% filter(scaffold %in% scaffs_12)

syntenyplot <- plot_chr12(ssr_12, window_width = 2e5+1)

ggsave("fig_chr12_synteny.pdf",syntenyplot, width=8, height=5, useDingbats=FALSE)


ggplot(total_win_df)+ geom_point(aes(x=pi_716, y=sites))#+coord_cartesian(xlim=c(0,0.01),ylim=c(0, 0.01))

#####
## kmer plots
#####
top_kmer_plot <- ggplot()+
  geom_line(data=filter(top_rand_kmer,cat=='random'), aes(x=x, y=kpct*100, group =samp), alpha =0.5, color="grey80")+
  geom_line(data=filter(top_rand_kmer,cat!='random'), aes(x=x, y=kpct*100, color =cat), size=1.2)+
  geom_point(data=filter(top_rand_kmer,cat!='random'), aes(x=x, y=kpct*100, color =cat), size=2)+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = c(0.85, 0.85), legend.background = element_blank())+
  labs(x="10Kb windows ranked by sex-specific kmer content", y="% of k-mers mapped to window")
d_kmer <- rand_kmer %>%
  mutate(nk = map_int(windowed, ~length(.$kpct)))%>%
  mutate(unmapped = map_dbl(windowed, ~filter(.x, id =="*_1")%>%summarise(kpct=sum(kpct) )%>%pull(kpct) )) %>%
  select(samp, cat, unmapped, nk)%>%
  mutate(libsize = map_int(samp, get.totalsize))

ggplot(d_kmer)+
  geom_point(aes(y=nk, x= libsize, color=cat), alpha =0.5)+
  labs(x="Library size", y="Number of k-mers found")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'none')
ggplot(filter(d_kmer, libsize >4e8, libsize <6e8))+
  geom_point(aes(y=nk, x= libsize, color=cat), alpha =0.5)+
  labs(x="Library size", y="Number of k-mers found")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'none')


unmappdist_plot <- ggplot()+
  geom_histogram(data=filter(d_kmer,cat=='random'), aes(x=unmapped*100), alpha =0.5, bins=50)+
  geom_vline(data=filter(d_kmer,cat!='random'), aes(xintercept=unmapped*100, color =cat), size=1.5)+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'none')+
  labs(x="% of unmapped k-mers")

totalkdist_plot<- ggplot()+
  geom_histogram(data=filter(d_kmer,cat=='random'), aes(x=nk), alpha =0.5, bins=50)+
  geom_vline(data=filter(d_kmer,cat!='random'), aes(xintercept=nk, color =samp), size=1.5)+
  labs(x="Total k-mers found")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'none')

kmerdist_plot <-  (top_kmer_plot | (totalkdist_plot / unmappdist_plot))
ggsave("random_kmerdist.pdf", kmerdist_plot, width= 8, height=5)

#####
### Top scaffold stats
#####
# Fst and Dxy btw males and females in two pops
zoom_raw_MF716 <- read_csv(gzfile("FstDxy_cal/topscaff_716_sex-diff.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  rename(pi_M716 = "pi_popB", pi_F716 =  "pi_popA", dxy_MF716 = "dxy_popA_popB", fst_MF716 = "Fst_popA_popB" )

zoom_raw_MF670 <- read_csv(gzfile("FstDxy_cal/topscaff_670_sex-diff.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  transmute(id, pi_M670 = pi_popB, pi_F670 = pi_popA, dxy_MF670 = dxy_popA_popB, fst_MF670 = Fst_popA_popB)

# Fst and Dxy btw two pops
zoom_raw_popdiff <- read_csv(gzfile("FstDxy_cal/topscaff_670_vs_716.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  transmute(id, pi_670 = pi_popA, pi_716 =  pi_popB)

zoom_raw_m_popdiff <- read_csv(gzfile("FstDxy_cal/topscaff_m_pop_diff.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  rename(pi_670 = "pi_popA", pi_716 =  "pi_popB", popDxy = "dxy_popA_popB", popFst = "Fst_popA_popB" )
zoom_raw_f_popdiff <- read_csv(gzfile("FstDxy_cal/topscaff_f_pop_diff.HD5.csv.gz"), col_types = "ciiiidddd")%>%
  unite("id", scaffold, start, sep = "_", remove = F)%>%
  rename(pi_670 = "pi_popA", pi_716 =  "pi_popB", popDxy = "dxy_popA_popB", popFst = "Fst_popA_popB" )

zoom_popdiff <- zoom_raw_m_popdiff %>%  
  transmute(id, 
            popDxy = (popDxy+ zoom_raw_f_popdiff$popDxy)/2, 
            popFst=(popFst+ zoom_raw_f_popdiff$popFst)/2 )%>%
  full_join(zoom_raw_popdiff, by='id')

# merge the three dataframes
zoom_pgstats <- full_join(zoom_raw_MF716, zoom_raw_MF670, by='id' ) %>% left_join(zoom_popdiff)%>%
  mutate(popFst = if_else(popFst < 0, 0, popFst),
         fst_MF670 = if_else(fst_MF670 < 0, 0, fst_MF670),
         fst_MF716 = if_else(fst_MF716 < 0, 0, fst_MF716))%>%
  mutate(sexDxy = (dxy_MF670 + dxy_MF716)/2 )%>%
  mutate(sexFst = (fst_MF670 + fst_MF716)/2)%>%
  mutate( piF = 0.5*(pi_F670 + pi_F716),
          piM = 0.5*(pi_M670 + pi_M716),
          pipop =0.5*(pi_670 + pi_716) )

zoom_pgstats_win  <-  win_1k %>% filter(scaffold %in% top_export$scaffold)%>%
  left_join(zoom_pgstats, by=c('id', 'scaffold', 'start')) %>%
  mutate(end=start+1e3)%>%
  replace_na(list(sexFst=0, piF=0, piM=0, sexDxy=0))

zoom_total_win_df <- full_join(zoom_pgstats_win, ssr_clean, by=c("id", "scaffold", "start")) 

