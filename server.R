library("shiny")
library("genefilter")
library("gplots")
library("DESeq2")
library("BiocParallel")

shinyServer(function(input, output) {
    
    # Access uploaded file if its supplied
    # .. else read in downloaded data
    dataTable <- reactive({
        inFile <- input$file1
        if (!is.null(inFile)) {
            myData <- read.table(inFile$datapath,header=T,sep="\t",row.names=1)
        } else {
            URL = "http://bowtie-bio.sourceforge.net/recount/countTables/sultan_count_table.txt"
            myData <- read.table(URL, header=T, sep="\t", nrows=1000, row.names=1)
        }
        return(myData)
    })
    
    # User Selection (2 x sliders & 2 x selectize menu)
    output$P_Value <- renderUI({
        myData <- dataTable()
        sliderInput("decimal",
                    "P-Value significance:",
                    min = 0,
                    max = 0.005,
                    step = 0.0005,
                    value = 0.002)
    })
    output$FoldChange <- renderUI({
        myData <- dataTable()
        sliderInput("integer",
                    "Differential gene fold change:",
                    min = 0,
                    max = 20,
                    step = 0.5,
                    value = 5)
    })
    output$TopGenes <- renderUI({
        myData <- dataTable()
        sliderInput("top",
                    "Select number of genes to visualise:",
                    min = 0,
                    max = 500,
                    step = 10,
                    value = 200)
    })
    output$TestSet <- renderUI({ 
        myData <- dataTable()
        selectInput(inputId="test",
                    label="Test Set:",
                    choices=names(myData), 
                    selected=names(myData)[1:2],
                    multiple=TRUE, 
                    selectize=TRUE)  
    })
    output$ControlSet <- renderUI({ 
        myData <- dataTable()
        selectInput(inputId="ctrl",
                    label="Control Set:",
                    choices=names(myData), 
                    selected=names(myData)[3:4],
                    multiple=TRUE, 
                    selectize=TRUE)
    })
    
    # Heatmap Plot
    output$Heatmap <- renderPlot({
        myData <- dataTable()
        logMatrix <- as.matrix(log(myData+1))
        sampleTable <- data.frame(samples=names(myData), condition=c("test","test","ctrl","ctrl"))
        colnames(sampleTable) <- c("sample","condition")
        heatCols <- colorRampPalette(c("blue","grey","red"))(n=100)
        sideCols <- c("darkgreen","black")[ sampleTable$condition ]
        rowvars <- rowVars(logMatrix)
        select = order(-rowvars)[1:input$top]
        logMatrix_top = logMatrix[select,]
        logMatrix_top_norm = t(scale(t(logMatrix_top),center = T, scale = T))
        heatmap.2(logMatrix_top_norm, margins=c(12,8),
                  key=TRUE, symkey=FALSE, 
                  Colv = T, Rowv = T, labRow=F,
                  ColSideColors=sideCols,
                  col = heatCols, density.info="none",
                  trace="none", scale="none",
                  dendrogram = "both", srtCol=45,
                  main=paste("Heatmap of top ",input$top," variable genes")
        )
    })

    # Data Contents - limited to 10 rows per tabset
    output$content <- renderDataTable({
        myData <- dataTable()
        },options = list(orderClasses = TRUE, 
                     autoWidth=TRUE, 
                     pageLength = 10)
    )
    
    # Differential Gene Expression using DESeq2
    dataInput <- reactive({
        myData <- dataTable()
        sampleTable <- data.frame(samples=names(myData), 
                                  condition=c("test","test","ctrl","ctrl"))
        colnames(sampleTable) <- c("sample","condition")
        dds <- DESeqDataSetFromMatrix(countData = myData, 
                                      colData = sampleTable,
                                      design = ~ condition)
        rm <- rowMeans(counts(dds))
        dds <- estimateSizeFactors(dds,controlGenes=which(rm>100))
        register(MulticoreParam(2)) # set this to however many cores are available
        vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
        vsdMatrix <- assay(vsd)
        cpmNorm <- counts(dds, normalized=T)
        colnames(cpmNorm) <- dds$condition
        dds <- DESeq(dds, fitType='local', parallel = TRUE)
        
        results <- results(dds, contrast=c("condition","test","ctrl"))
        results_ord <- results[order(results$padj),]
    })
    
    output$MA_plot <- renderPlot({
        plotMA(dataInput(), alpha= input$decimal, 
               ylim=c(-input$integer,input$integer),
               main=paste("MA Plot of Differential Gene Expression (logFC=",
                          input$integer,", p-value=",input$decimal,")")
               )
    })
}) 