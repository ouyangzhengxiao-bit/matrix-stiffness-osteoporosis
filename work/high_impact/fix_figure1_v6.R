suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(scales)})

root <- normalizePath('.')
out <- file.path(root, 'outputs/研究升级/Communications_Biology_v6')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_classic(base_size = 9, base_family = 'Arial'))
cols <- c(blue='#2C6EAA', orange='#D87A2C', grey='#9AA1A6')

# a: analysis architecture. Use explicit box boundaries so connector arrows
# terminate in the whitespace between nodes rather than entering the boxes.
nodes <- data.frame(
  x = c(.90, 2.25, 3.68, 5.12, 6.55), y = 1,
  half_w = c(.40, .45, .49, .40, .49),
  half_h = c(.25, .21, .25, .21, .21),
  label = c('7 independent\nGEO studies\n13 contrasts',
            'Study-clustered\ngene effects',
            '12,241 genes\nrandom-effects\nmeta-analysis',
            'Human bone\ncell programs',
            'ATAC and eBMD\nfalsification'))
arrows <- data.frame(
  x = nodes$x[-nrow(nodes)] + nodes$half_w[-nrow(nodes)] + .08,
  xend = nodes$x[-1] - nodes$half_w[-1] - .12,
  y = 1, yend = 1)
p1a <- ggplot() +
  geom_segment(data=arrows, aes(x=x,xend=xend,y=y,yend=yend),
               arrow=arrow(length=unit(1.7,'mm'),type='closed'),
               linewidth=.5, colour='#555555') +
  geom_rect(data=nodes,
            aes(xmin=x-half_w,xmax=x+half_w,ymin=y-half_h,ymax=y+half_h),
            fill='white',colour='black',linewidth=.35) +
  geom_text(data=nodes,aes(x=x,y=y,label=label),size=2.55,lineheight=.94) +
  annotate('text', x=.58, y=1.42, label='a', fontface='bold', size=4) +
  coord_cartesian(xlim=c(.42,7.12), ylim=c(.68,1.42), clip='off') +
  theme_void()

# b: study-level gene-effect concordance
cor <- read.csv(file.path(root,'outputs/研究升级/meta_analysis_v2/study_effect_spearman.csv'), check.names=FALSE)
rownames(cor) <- cor[,1]; cor <- cor[,-1]
hm <- as.data.frame(as.table(as.matrix(cor))); names(hm) <- c('Study1','Study2','rho')
p1b <- ggplot(hm, aes(Study1,Study2,fill=rho)) +
  geom_tile(colour='white', linewidth=.25) +
  scale_fill_gradient2(low=cols[['blue']], mid='white', high=cols[['orange']],
                       limits=c(-.26,.26), name='Spearman\nrho') +
  coord_equal() +
  labs(x=NULL,y=NULL,title='b  Gene-effect concordance') +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6.5),
        axis.text.y=element_text(size=6.5),
        plot.title=element_text(face='bold',size=9),
        plot.margin=margin(3,4,2,2))

# c: inference-unit sensitivity. Wrap long pathway names and allocate more width.
pooled <- read.csv(file.path(root,'outputs/研究升级/meta_analysis_v2/pathway_meta_analysis.csv'))
study <- read.csv(file.path(root,'outputs/研究升级/meta_analysis_v2/pathway_study_unit_meta.csv'))
top <- pooled[order(pooled$BH_q),][1:6,]
z <- merge(top[,c('pathway','BH_q')], study[,c('pathway','BH_q_study_unit')], by='pathway')
z$pathway <- factor(z$pathway, levels=rev(z$pathway))
zz <- rbind(
  data.frame(pathway=z$pathway, method='Pooled gene ranks', q=z$BH_q),
  data.frame(pathway=z$pathway, method='Independent studies', q=z$BH_q_study_unit))
zz$label <- gsub('REACTOME_|HALLMARK_','',zz$pathway)
zz$label <- gsub('_',' ',zz$label)
zz$label <- vapply(zz$label, function(x) paste(strwrap(x, width=31), collapse='\n'), character(1))
zz$label <- factor(zz$label, levels=rev(unique(as.character(zz$label))))
p1c <- ggplot(zz, aes(-log10(pmax(q,1e-14)),label,colour=method,shape=method)) +
  geom_vline(xintercept=-log10(.05),linetype=2,colour='#777777') +
  geom_point(size=2.3) +
  scale_colour_manual(values=c('Pooled gene ranks'=cols[['orange']],
                               'Independent studies'=cols[['blue']])) +
  labs(x=expression(-log[10]('FDR q')), y=NULL,
       title='c  Inference-unit sensitivity', colour=NULL, shape=NULL) +
  guides(colour=guide_legend(nrow=1,byrow=TRUE), shape=guide_legend(nrow=1,byrow=TRUE)) +
  theme(axis.text.y=element_text(size=6.3,lineheight=.9),
        legend.position='bottom', legend.direction='horizontal',
        legend.text=element_text(size=7), legend.key.width=unit(4,'mm'),
        legend.spacing.x=unit(2,'mm'),
        plot.title=element_text(face='bold',size=9),
        plot.margin=margin(3,8,0,2))

fig1 <- p1a / (p1b | p1c) +
  plot_layout(heights=c(.31,1), widths=c(.86,1.30)) &
  theme(plot.margin=margin(3,7,2,3))

ggsave(file.path(out,'Figure_1.png'), fig1, width=8.2, height=4.75, dpi=600, bg='white')
ggsave(file.path(out,'Figure_1.tif'), fig1, width=8.2, height=4.75, dpi=600,
       compression='lzw', bg='white')
ggsave(file.path(out,'Figure_1.pdf'), fig1, width=8.2, height=4.75,
       device=cairo_pdf, bg='white')
cat('Corrected Figure 1 written to', out, '\n')
