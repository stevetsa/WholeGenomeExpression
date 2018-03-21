library(shiny)
library(shinyjs)
library(DT)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(reshape2)
library(plotly)
library(Rsubread)

#setwd('/Users/valdezkm/Documents/Hackathon')
load('./viSRA_practice_data.RData')

system(paste0('/home/ubuntu/ncbi-magicblast-1.3.0/bin/makeblastdb -in /home/ubuntu/hgDir/chr21.fa -dbtype nucl -parse_seqids -out Chr21 -title "Chr21"'))


shinyServer(function(input, output) {
  observeEvent(
    input$SRAbutton,
    isolate({ 
      write.csv(expression, file='normalized_data.csv')
      expression
    })
  )
  observeEvent(
    input$SRAbutton,
    isolate({
      SRA1 = input$SRRcode1
      SRA2 = input$SRRcode2
      
      #### Blast SRR numbers against reference genome ####
      #system(paste0(magicDir,'/magicblast -sra ',SRA1,SRA2,' -db ',magicDir,'/GRCh38'))
      
      #### Count reads with Rsubread ####
      # subreadDir = '/home/ubuntu/subread-1.6.0-Linux-x86_64/bin'
      # GTF = 'gencode.v21.annotation.gtf'
      # 
      # counts1 = system(paste0(subreadDir,'/featureCounts -t exon -g gene_id -a ',GTF,' -o counts1.txt SRA1.bam'))
      # counts2 = system(paste0(subreadDir,'/featureCounts -t exon -g gene_id -a ',GTF,' -o counts2.txt SRA2.bam'))
      # 

      
      #write.table(SRA1, file = 'sra_cond1', sep = "\n", row.names = F, quote = F, col.names = F)
      #write.table(SRA1, file = 'sra_cond2', sep = "\n", row.names = F, quote = F, col.names = F)
    })
  )
  observeEvent(input$SRAbutton,
    isolate({
    dat <- as.data.frame(expression)
    listGenes = read.delim(input$listOfGenes$datapath, sep = '\n', header = F)
    for (i in 1:length(listGenes[,1])) {
      gene = listGenes[i,1]
      jpeg(file = paste0(gene,".jpeg"))
      print(dat %>%
        mutate(geneID = rownames(dat)) %>%
        filter(geneID == toupper(gene)) %>%
        melt() %>%
        ggplot(aes(x = variable, y = value)) +
        geom_point() +
        theme_bw() +
        xlab("SRA ID") +
        ylab("Expression (TPM)") +
        ggtitle(paste(gene)))
      dev.off()
    }
    })
  )
  output$sra=DT::renderDataTable(DT::datatable(
    {
    expression
    })
  )
  output$ssgsea=DT::renderDataTable(DT::datatable(
    {
      geneSet =  getGmt(input$geneSet)
      ssResults = gsva(as.matrix(expression),geneSet,method='ssgsea')
      write.csv(ssResults, file = 'ssGSEA_pathways.csv')
    })
  )
  output$dotPlot <- renderPlotly({
    input$dotPlotButton
    isolate({
      dat <- as.data.frame(expression)
      d <- dat %>%
        mutate(geneID = rownames(dat)) %>%
        filter(geneID == toupper(input$geneName)) %>%
        melt() %>%
        ggplot(aes(x = variable, y = value)) +
        geom_point() +
        theme_bw() +
        xlab("SRA ID") +
        ylab("Expression (TPM)") +
        ggtitle(paste(input$geneName))
      d <- plotly_build(d)
      d$elementId <- NULL
      print(d)
    })
  })
  output$inc = renderUI(
    {
      tags$iframe(
        src='testing_fastqc.html',
        width='100%',
        height='800px')
    }
  )
})

