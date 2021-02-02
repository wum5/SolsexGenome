#!/Users/mengwu/Documents/Research/CAFE_v3/cafe/cafe
#version
#date
tree (((((lycopersicum:8,tuberosum:8):2,appendiculatum:10):7,sinuosa:17):7,attenuata:24):6,(inflata:3,axillaris:3):27)
load -i filtered_cafe_input.txt -p 0.01 -t 2 -l log.txt -filter
lambda -s
report resultfile

