#### PROJECT PROGRESSION MODEL USING TRONCO PACKAGE  


library(devtools)
library(pheatmap)
library(gridExtra)
library(vioplot)
library(xlsx)
library(energy)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)


library(TCGAbiolinks)
library(TRONCO)

#### INPUT DATA ####
# Import MAF data set 
maf <- GDCquery_Maf("SKCM", pipelines = "mutect2")
# save as a data-frame
maf<- as.data.frame(maf)

## analyze the data using TRONCO package 
# import the maf file into the TRONCO package 
maf_data = import.MAF(maf, merge.mutation.types = T)
view(maf_data) # make sure the data was converted properly 

maf_data<- events.selection(as.alterations(maf_data), filter.freq = .05)

# filter the data to hypothesized events based on literature review/ cbioportal etc.  
drivers<- c('BRAF', 'NRAS', 'RASA2', 'CTLA4', 'MUC20',  'PPP6C', 'RAC1', 'PTEN', 'KMT2C', 'TP53', 'CKIT','TP63', 'NF1','BCL2' , 'KDR', 'CDKN2A','AGO2', 'KMT2A', 'RAC1', 'ARID2', 'TERT', 'PI3K', 'KIT', 'MET', 'MAP2K1', 'MAPK', 'C288A' )

# exclusivity events
gene_hypothesis<- hypothesis.add(gene_hypothesis, 'NRAS xor BRAF', XOR('NRAS', 'BRAF'))

skcm_dataset<- events.selection(maf_data, filter.in.names = drivers)

# change colors, dark green is too dark, not sure why events wont further seperated???
skcm_dataset<- change.color(skcm_dataset, 'Mutation', 'purple')

#print the oncoprint
oncoprint(skcm_dataset, title = 'Select Melanoma Mutations',
          samples.cluster = T,
          genes.cluster = T,
          gene.annot.color = 'Set2'
          )

#### Run the Capri Model ####
model<-tronco.capri(skcm_dataset,
             nboot = 5)

## plot the model 
tronco.plot(model, 
            fontsize = 14, 
            scale.nodes = .45, 
            height.logic = 0.5, 
            confidence = 'tp',
            legend.cex = 0.7, 
            legend.pos = 'top',
            label.edge.size = 12,
            title = 'Melanoma CAPRI Model')

view(model)


