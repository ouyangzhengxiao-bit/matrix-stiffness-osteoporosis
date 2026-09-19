suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
root <- normalizePath('.')
out <- file.path(root,'outputs/研究升级/Communications_Biology_v6')
theme_set(theme_classic(base_size=8,base_family='Arial'))
cols <- c(orange='#D87A2C',teal='#25877A',grey='#9AA1A6',purple='#7B61A8')

ind <- read.csv(file.path(root,'outputs/研究升级/high_impact_v4/independent_osteoblast_validation/independent_osteoblast_study_effects.csv'))
ind$study <- factor(ind$study,levels=ind$study[order(ind$raw_matched_effect)])
p3a <- ggplot(ind,aes(raw_matched_effect,study))+
  geom_vline(xintercept=0,colour='#777777',linewidth=.4)+
  geom_point(size=2.5,colour=cols[['orange']])+
  labs(x='Matched difference in mean absolute study effect',y=NULL,
       title='a  Independent human osteoblast program',
       subtitle='GSE147390 markers: 6/7 studies positive; exact sign-flip P=0.039')+
  theme(plot.title=element_text(face='bold',size=8),plot.subtitle=element_text(size=7))

atac_json <- fromJSON(file.path(root,'outputs/研究升级/high_impact_v4/atac_validation/posthoc_tpm_sensitivity_summary.json'))
atac <- data.frame(comparison=factor(c('Mechanical stretch\nForce+ vs senescent','Soft matrix\nForce- vs young'),levels=c('Mechanical stretch\nForce+ vs senescent','Soft matrix\nForce- vs young')),
 global=c(atac_json$force_plus$genomewide_mean_log1p_TPM_delta,atac_json$soft_force_minus$genomewide_mean_log1p_TPM_delta),
 specific=c(atac_json$force_plus$directional_score,atac_json$soft_force_minus$directional_score))
aa <- reshape(atac,direction='long',varying=c('global','specific'),v.names='effect',timevar='endpoint',times=c('Genome-wide promoter shift','Signature directional contrast'))
p3b <- ggplot(aa,aes(comparison,effect,fill=endpoint))+
  geom_hline(yintercept=0,colour='#777777')+geom_col(position=position_dodge(width=.7),width=.62)+
  scale_fill_manual(values=c('Genome-wide promoter shift'=cols[['teal']],'Signature directional contrast'=cols[['purple']]))+
  labs(x=NULL,y='Mean log1p TPM change',title='b  ATAC response',subtitle='No signature specificity',fill=NULL)+
  theme(axis.text.x=element_text(size=6.5),legend.position='bottom',
        plot.title=element_text(face='bold',size=8),plot.subtitle=element_text(size=7),
        legend.text=element_text(size=6.5))

gen_json <- fromJSON(file.path(root,'outputs/研究升级/high_impact_v4/genetics_validation/genetics_validation_summary.json'))
ov <- gen_json$significant_eBMD_overlap
over <- data.frame(category=c('Observed','Expected under overlap null'),genes=c(ov$overlap,ov$signature*ov$eBMD_significant/ov$universe))
p3c <- ggplot(over,aes(category,genes,fill=category))+
  geom_col(width=.55)+geom_text(aes(label=sprintf('%.1f',genes)),vjust=-.4,size=2.7)+
  scale_fill_manual(values=c('Observed'=cols[['orange']],'Expected under overlap null'=cols[['grey']]))+
  coord_cartesian(ylim=c(0,40))+
  labs(x=NULL,y='eBMD-significant genes',title='c  eBMD overlap',
       subtitle=sprintf('No enrichment; P=%.3f',ov$hypergeometric_p))+
  theme(legend.position='none',axis.text.x=element_text(size=6.5),
        plot.title=element_text(face='bold',size=8),plot.subtitle=element_text(size=7))

fig3 <- p3a/(p3b|p3c)+plot_layout(heights=c(1.18,1))
ggsave(file.path(out,'Figure_3.png'),fig3,width=7.2,height=6.2,dpi=600,bg='white')
ggsave(file.path(out,'Figure_3.tif'),fig3,width=7.2,height=6.2,dpi=600,compression='lzw',bg='white')
ggsave(file.path(out,'Figure_3.pdf'),fig3,width=7.2,height=6.2,device=cairo_pdf,bg='white')
cat('Corrected Figure 3 written to',out,'\n')
