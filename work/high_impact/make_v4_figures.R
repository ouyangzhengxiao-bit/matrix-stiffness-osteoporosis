suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(scales); library(jsonlite)})

root <- normalizePath('.')
out <- file.path(root,'outputs/研究升级/Communications_Biology_v4')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
theme_set(theme_classic(base_size=9, base_family='Arial'))
cols <- c(blue='#2C6EAA', orange='#D87A2C', teal='#25877A', grey='#9AA1A6', red='#B04A4A', purple='#7B61A8')

# Figure 1: architecture, discordance, and inference-unit contrast
nodes <- data.frame(x=c(1,2.3,3.6,5.0,6.4), y=1,
 label=c('7 independent\nGEO studies\n13 contrasts','Study-clustered\ngene effects','12,241 genes\nrandom-effects\nmeta-analysis','Human bone\ncell programs','ATAC and eBMD\nfalsification'))
p1a <- ggplot(nodes,aes(x,y)) +
 geom_segment(data=data.frame(x=c(1.35,2.65,3.95,5.35),xend=c(1.95,3.25,4.65,6.05),y=1,yend=1),
              aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(length=unit(2.2,'mm')),linewidth=.5,colour='#555555')+
 geom_label(aes(label=label),size=2.6,linewidth=.35,label.padding=unit(.18,'lines'),fill='white')+
 annotate('text',x=.6,y=1.55,label='a',fontface='bold',size=4)+
 coord_cartesian(xlim=c(.5,6.9),ylim=c(.55,1.55),clip='off')+theme_void()

cor <- read.csv(file.path(root,'outputs/研究升级/meta_analysis_v2/study_effect_spearman.csv'),check.names=FALSE)
rownames(cor)<-cor[,1]; cor<-cor[,-1]
hm <- as.data.frame(as.table(as.matrix(cor))); names(hm)<-c('Study1','Study2','rho')
p1b <- ggplot(hm,aes(Study1,Study2,fill=rho))+geom_tile(colour='white',linewidth=.25)+
 scale_fill_gradient2(low=cols[['blue']],mid='white',high=cols[['orange']],limits=c(-.26,.26),name='Spearman\nrho')+
 coord_equal()+labs(x=NULL,y=NULL,title='b  Gene-effect concordance')+
 theme(axis.text.x=element_text(angle=45,hjust=1,size=6.5),axis.text.y=element_text(size=6.5),plot.title=element_text(face='bold',size=9))

pooled <- read.csv(file.path(root,'outputs/研究升级/meta_analysis_v2/pathway_meta_analysis.csv'))
study <- read.csv(file.path(root,'outputs/研究升级/meta_analysis_v2/pathway_study_unit_meta.csv'))
top <- pooled[order(pooled$BH_q),][1:6,]
z <- merge(top[,c('pathway','BH_q')],study[,c('pathway','BH_q_study_unit')],by='pathway')
z$pathway <- factor(z$pathway,levels=rev(z$pathway))
zz <- rbind(data.frame(pathway=z$pathway,method='Pooled gene ranks',q=z$BH_q),
            data.frame(pathway=z$pathway,method='Independent studies',q=z$BH_q_study_unit))
zz$label <- gsub('REACTOME_|HALLMARK_','',zz$pathway); zz$label<-gsub('_',' ',zz$label)
zz$label<-factor(zz$label,levels=rev(unique(as.character(zz$label))))
p1c <- ggplot(zz,aes(-log10(pmax(q,1e-14)),label,colour=method,shape=method))+
 geom_vline(xintercept=-log10(.05),linetype=2,colour='#777777')+geom_point(size=2.3)+
 scale_colour_manual(values=c('Pooled gene ranks'=cols[['orange']],'Independent studies'=cols[['blue']]))+
 labs(x=expression(-log[10]('FDR q')),y=NULL,title='c  Inference-unit sensitivity',colour=NULL,shape=NULL)+
 theme(axis.text.y=element_text(size=6.5),legend.position='bottom',plot.title=element_text(face='bold',size=9))
fig1 <- p1a / (p1b | p1c) + plot_layout(heights=c(.38,1))

# Figure 2: human bone localization and study-unit audit
sc <- read.csv(file.path(root,'outputs/研究升级/high_impact_v4/human_bone_localization/human_scRNA_cluster_tests.csv'))
sc <- subset(sc,endpoint=='abs_meta_z'); sc$lineage <- ifelse(grepl('MSC|osteoblast',sc$cluster,ignore.case=TRUE),'Osteoblast lineage','Other')
sc$short <- gsub('^[0-9]+\\. ','',sc$cluster); sc$short<-gsub('Mature neutrophils \\(([^)]+)\\)','Neutrophils \\(\\1\\)',sc$short)
sc$short<-factor(sc$short,levels=sc$short[order(sc$effect)])
p2a <- ggplot(sc,aes(effect,short,colour=lineage))+geom_vline(xintercept=0,colour='#888888',linewidth=.35)+geom_point(size=2)+
 scale_colour_manual(values=c('Osteoblast lineage'=cols[['orange']],'Other'=cols[['grey']]))+
 labs(x='Matched difference in mean |meta z|',y=NULL,title='a  Human femoral single-cell programs',colour=NULL)+
 theme(axis.text.y=element_text(size=6.2),legend.position='bottom',plot.title=element_text(face='bold'))

sp <- read.csv(file.path(root,'outputs/研究升级/high_impact_v4/human_bone_localization/human_spatial_cluster_tests.csv'))
sp <- subset(sp,endpoint=='abs_meta_z'); sp$lineage<-ifelse(grepl('MSC|Osteoblast',sp$cluster,ignore.case=TRUE),'Osteoblast lineage','Other')
sp$short<-sp$cluster; sp$short<-factor(sp$short,levels=sp$short[order(sp$effect)])
p2b <- ggplot(sp,aes(effect,short,colour=lineage))+geom_vline(xintercept=0,colour='#888888',linewidth=.35)+geom_point(size=2)+
 scale_colour_manual(values=c('Osteoblast lineage'=cols[['orange']],'Other'=cols[['grey']]))+
 labs(x='Matched difference in mean |meta z|',y=NULL,title='b  Human bone spatial programs',colour=NULL)+
 theme(axis.text.y=element_text(size=6.5),legend.position='bottom',plot.title=element_text(face='bold'))

su <- read.csv(file.path(root,'outputs/研究升级/high_impact_v4/human_bone_localization/study_unit_localization_effects.csv'))
su$program_source<-factor(su$program_source,levels=c('scRNA','spatial'),labels=c('Single-cell program','Spatial program'))
p2c <- ggplot(su,aes(raw_matched_effect,reorder(study,raw_matched_effect),colour=program_source))+
 geom_vline(xintercept=0,colour='#777777',linewidth=.35)+geom_point(size=2)+facet_wrap(~program_source,ncol=2,scales='free_y')+
 scale_colour_manual(values=c('Single-cell program'=cols[['orange']],'Spatial program'=cols[['teal']]))+
 labs(x='Study-level matched enrichment effect',y=NULL,title='c  Independent-study robustness')+
 theme(legend.position='none',plot.title=element_text(face='bold'))
fig2 <- (p2a | p2b) / p2c + plot_layout(heights=c(1.25,.85))

# Figure 3: multi-layer tests of transferability
ind <- read.csv(file.path(root, "outputs/研究升级/high_impact_v4/independent_osteoblast_validation/independent_osteoblast_study_effects.csv"))
ind$study <- factor(ind$study, levels = ind$study[order(ind$raw_matched_effect)])
p3a <- ggplot(ind, aes(raw_matched_effect, study)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.4) +
  geom_point(size = 2.7, colour = cols[["orange"]]) +
  labs(x = "Matched difference in mean absolute study effect", y = NULL,
       title = "a  Independent human osteoblast program",
       subtitle = "GSE147390 markers: 6/7 studies positive; exact sign-flip P=0.039") +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8))

atac_json <- fromJSON(file.path(root,'outputs/研究升级/high_impact_v4/atac_validation/posthoc_tpm_sensitivity_summary.json'))
atac <- data.frame(comparison=factor(c('Mechanical stretch\nForce+ vs senescent','Soft matrix\nForce- vs young'),levels=c('Mechanical stretch\nForce+ vs senescent','Soft matrix\nForce- vs young')),
 global=c(atac_json$force_plus$genomewide_mean_log1p_TPM_delta,atac_json$soft_force_minus$genomewide_mean_log1p_TPM_delta),
 specific=c(atac_json$force_plus$directional_score,atac_json$soft_force_minus$directional_score))
aa <- reshape(atac,direction='long',varying=c('global','specific'),v.names='effect',timevar='endpoint',times=c('Genome-wide promoter shift','Signature directional contrast'))
p3b<-ggplot(aa,aes(comparison,effect,fill=endpoint))+geom_hline(yintercept=0,colour='#777777')+geom_col(position=position_dodge(width=.7),width=.62)+
 scale_fill_manual(values=c('Genome-wide promoter shift'=cols[['teal']],'Signature directional contrast'=cols[['purple']]))+
 labs(x=NULL,y='Mean log1p TPM change',title='b  Global ATAC shift without signature specificity',fill=NULL)+
 theme(axis.text.x=element_text(size=7),legend.position='bottom',plot.title=element_text(face='bold'))

gen_json <- fromJSON(file.path(root,'outputs/研究升级/high_impact_v4/genetics_validation/genetics_validation_summary.json'))
ov <- gen_json$significant_eBMD_overlap
over <- data.frame(category=c('Observed','Expected under overlap null'),genes=c(ov$overlap,ov$signature*ov$eBMD_significant/ov$universe))
p3c<-ggplot(over,aes(category,genes,fill=category))+geom_col(width=.55)+geom_text(aes(label=sprintf('%.1f',genes)),vjust=-.4,size=3)+
 scale_fill_manual(values=c('Observed'=cols[['orange']],'Expected under overlap null'=cols[['grey']]))+coord_cartesian(ylim=c(0,40))+
 labs(x=NULL,y='eBMD-significant genes',title=sprintf('c  No eBMD enrichment\nP=%.3f',ov$hypergeometric_p))+
 theme(legend.position='none',axis.text.x=element_text(size=7),plot.title=element_text(face='bold',size=9))
fig3 <- p3a / (p3b | p3c) + plot_layout(heights=c(1.2,1))

save_all <- function(plot,name,w=7.2,h=6.0){
 ggsave(file.path(out,paste0(name,'.png')),plot,width=w,height=h,dpi=320,bg='white')
 ggsave(file.path(out,paste0(name,'.tif')),plot,width=w,height=h,dpi=320,compression='lzw',bg='white')
 ggsave(file.path(out,paste0(name,'.pdf')),plot,width=w,height=h,device=cairo_pdf,bg='white')
}
save_all(fig1,'Figure_1',7.2,5.7); save_all(fig2,'Figure_2',7.2,7.0); save_all(fig3,'Figure_3',7.2,6.2)
cat('Figures written to',out,'\n')
