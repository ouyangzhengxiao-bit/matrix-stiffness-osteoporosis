suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

root <- normalizePath('.')
out <- file.path(root, 'outputs/研究升级/Communications_Biology_v6')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
theme_set(theme_classic(base_size=9, base_family='Arial'))
cols <- c(orange='#D87A2C', teal='#25877A', grey='#9AA1A6')

sc <- read.csv(file.path(root,'outputs/研究升级/high_impact_v4/human_bone_localization/human_scRNA_cluster_tests.csv'))
sc <- subset(sc, endpoint=='abs_meta_z')
sc$lineage <- ifelse(grepl('MSC|osteoblast',sc$cluster,ignore.case=TRUE),'Osteoblast lineage','Other')
sc$short <- gsub('^[0-9]+\\. ','',sc$cluster)
sc$short <- gsub('Mature neutrophils \\(([^)]+)\\)','Neutrophils \\(\\1\\)',sc$short)
sc$short <- factor(sc$short, levels=sc$short[order(sc$effect)])
p2a <- ggplot(sc,aes(effect,short,colour=lineage)) +
  geom_vline(xintercept=0,colour='#888888',linewidth=.35) + geom_point(size=2) +
  scale_colour_manual(values=c('Osteoblast lineage'=cols[['orange']],'Other'=cols[['grey']])) +
  labs(x='Matched difference in mean |meta z|',y=NULL,
       title='a  Human femoral single-cell programs',colour=NULL) +
  guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
  theme(axis.text.y=element_text(size=6.2),legend.position='bottom',
        legend.text=element_text(size=7),plot.title=element_text(face='bold',size=9),
        plot.margin=margin(3,5,2,3))

sp <- read.csv(file.path(root,'outputs/研究升级/high_impact_v4/human_bone_localization/human_spatial_cluster_tests.csv'))
sp <- subset(sp,endpoint=='abs_meta_z')
sp$lineage <- ifelse(grepl('MSC|Osteoblast',sp$cluster,ignore.case=TRUE),'Osteoblast lineage','Other')
sp$short <- factor(sp$cluster,levels=sp$cluster[order(sp$effect)])
p2b <- ggplot(sp,aes(effect,short,colour=lineage)) +
  geom_vline(xintercept=0,colour='#888888',linewidth=.35) + geom_point(size=2) +
  scale_colour_manual(values=c('Osteoblast lineage'=cols[['orange']],'Other'=cols[['grey']])) +
  labs(x='Matched difference in mean |meta z|',y=NULL,
       title='b  Human bone spatial programs',colour=NULL) +
  guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
  theme(axis.text.y=element_text(size=6.4),legend.position='bottom',
        legend.text=element_text(size=7),plot.title=element_text(face='bold',size=9),
        plot.margin=margin(3,10,2,4))

su <- read.csv(file.path(root,'outputs/研究升级/high_impact_v4/human_bone_localization/study_unit_localization_effects.csv'))
su$program_source <- factor(su$program_source,levels=c('scRNA','spatial'),
                            labels=c('Single-cell program','Spatial program'))
p2c <- ggplot(su,aes(raw_matched_effect,reorder(study,raw_matched_effect),colour=program_source)) +
  geom_vline(xintercept=0,colour='#777777',linewidth=.35) + geom_point(size=2) +
  facet_wrap(~program_source,ncol=2,scales='free_y') +
  scale_colour_manual(values=c('Single-cell program'=cols[['orange']],
                               'Spatial program'=cols[['teal']])) +
  labs(x='Study-level matched enrichment effect',y=NULL,
       title='c  Independent-study robustness') +
  theme(legend.position='none',plot.title=element_text(face='bold',size=9),
        strip.text=element_text(size=7.5),axis.text.y=element_text(size=6.5),
        plot.margin=margin(3,9,2,3))

fig2 <- (p2a | p2b) / p2c +
  plot_layout(heights=c(1.18,.82), widths=c(1.12,1)) &
  theme(plot.margin=margin(3,8,2,3))

ggsave(file.path(out,'Figure_2.png'),fig2,width=8.25,height=6.65,dpi=600,bg='white')
ggsave(file.path(out,'Figure_2.tif'),fig2,width=8.25,height=6.65,dpi=600,
       compression='lzw',bg='white')
ggsave(file.path(out,'Figure_2.pdf'),fig2,width=8.25,height=6.65,
       device=cairo_pdf,bg='white')
cat('Corrected Figure 2 written to',out,'\n')
