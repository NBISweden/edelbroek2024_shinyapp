# Shiny app for plotting RNA and protein data for Dictyostelium discoideum
#
# Jakub Orzechowski Westholm

library(shiny)
library(DT)
library(grid)
library(gridExtra)
library(tidyverse)

# Creates a scatterplot of gene expression vs. time, and fits a loess curve
# to the data 
# Input: 
# - The id of the gene to plot
# - A gene expression table
# - An (optional) array annotating the data, mapping row names in the gene expression table to strings.
# - A sample sheet, with column "Minute"
# - An unit to use in the y-axis label, e.g. "mRNA", "protein" etc.
# - A string to to use in the caption, e.g. "mRNA", "protein" etc.
# - Span parameter for fitting the loess curve (default 0.6)
plot.gene <- function(geneId, geneExp, geneImpute, sampleSheet, dataAnnot=NULL, yaxisUnit="", caption = "", span = 0.6, legendPosition="none") {
	sampleSheet <- sampleSheet[sampleSheet$Sample %in% colnames(geneExp),] # make sure samples match expression data

	plotData <- data.frame(expLevel = t(geneExp[geneId,]), 
												 imputed = t(geneImpute[geneId,]),
												 time = sampleSheet$Hours)
	colnames(plotData) <- c("expLevel", "imputed", "time")
	
	tmpSymbol <- ""
	if(!is.null(dataAnnot)){
		tmpSymbol <- dataAnnot[geneId]
		if(is.na(tmpSymbol)){tmpSymbol <- "" }
	}
	
	ggplot(data = plotData, aes(x = time, y = expLevel, color=imputed, label)) + geom_point( size=3) + 
		scale_color_manual(values=c("black", "gray70")) + 
		geom_smooth(method = "loess", span = span, fill="gray80", color="gray40") + 
		xlab("Time") + 
		ylab(yaxisUnit) +
		labs(title=caption,
				 subtitle=paste(geneId, tmpSymbol)) +
		theme_classic() +
		theme(axis.text.x = element_text(color="black", size = 14),
					axis.text.y = element_text(color="black", size = 14),
					axis.title.x = element_text(color="black", size=15),
					axis.title.y = element_text(color="black", size=15),
					plot.title = element_text(color="black", size=16, face="bold")) +
		theme(legend.position=legendPosition) +
		scale_x_continuous(
			breaks = c(0, 2, 4, 6, 8, 10),
			labels = c("0h", "2h", "4h", "6h", "8h", "10h")
		)
}


# Average gene expression over groups of samples (e.g. cell cycle phase)
# Input: 
# - A gene expression table
# - A list of arrays of column indexes, indicating which columns should be grouped together
# Output:
# - A new gene expression matrix, but with expression levels averged according to grouping given as input.
timePointAverages <- function(geneExp, timeGroups){
	sapply(timeGroups,
				 function(x){
				 	apply(geneExp[,intersect(x,colnames(geneExp))],1, mean, na.rm=T)
				 })
}



#############################
# Load data

#####
# Read gene annotations, downloaded from https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab
geneInfo <- read_tsv(file="shiny_data/gene_info.tsv", quote = "\"") %>%
	mutate(Gene_products = replace_na(Gene_products, "")) 
geneId <- geneInfo %>% as.data.frame() %>% .[,2]
names(geneId) <- geneInfo %>% as.data.frame() %>% .[,1]
geneSymbol <- geneInfo %>% as.data.frame() %>% .[,3]
names(geneSymbol) <- geneInfo %>% as.data.frame() %>% .[,1]
geneDescShort <- geneInfo %>% as.data.frame() %>% .[,4] 
names(geneDescShort) <- geneInfo %>% as.data.frame() %>% .[,1]


#####
# Load data
mRNAData <- read.table(file="shiny_data/rna_exp.tsv", sep="\t", header=T)
proteinData <- read.table(file="shiny_data/protein_exp.tsv", sep="\t", header=T)
proteinDataImputed <- read.table(file="shiny_data/protein_exp_is_imputed.tsv", sep="\t", header=T)
mRNADataInputed <- as.data.frame(mRNAData < 0) # data frame with only FALSE, since no RNA data is imputed

# Normalize so that the sum for each time point = the number of observations (RNAs or proteins)
norm_fun <- function(m){
	m <- as.matrix(m) 
	m_norm <- log2(t(t(2^m)/apply(2^m,2,sum)))
	return(as.data.frame(m_norm))
}
mRNAData <- norm_fun(mRNAData) 
proteinData <- norm_fun(proteinData) 


# Which samples to use, and which to leave out.
useRnaSamples <- colnames(mRNAData)
useProteinSamples <- colnames(proteinData)
useSamples <- sort(union(useRnaSamples, useProteinSamples))


sampleSheet <- data.frame(Sample=useSamples,
													Hours=sapply(colnames(mRNAData), function(x){as.numeric(gsub("h", "", strsplit(x, "_")[[1]][2]))}))
allGenes <- names(geneId)


#####
# Load and format results from differential expression analysis

# mRNA - EXPAND THIS TABLE?
limmaRna <- right_join(geneInfo, read_tsv(file="shiny_data/rna_de_res.tsv"), by="UniProt_ID") %>%
	dplyr::select(UniProt_ID, GENE_ID, Gene_Name, Gene_products, padj) %>%
	mutate(padj = signif(padj,2))  %>% 
	mutate(padj = replace_na(padj, 1)) %>% 
	arrange(padj) %>% 
	as.data.frame()
rownames(limmaRna) <- limmaRna$UniProt_ID

# proteomics - EXPAND THIS TABLE?
limmaProt <- right_join(geneInfo, read_tsv(file="shiny_data/protein_de_res.tsv"), by="UniProt_ID") %>%
	dplyr::select(UniProt_ID, GENE_ID, Gene_Name, Gene_products,adj.P.Val) %>%
	rename(padj = adj.P.Val) %>%
	mutate(padj = signif(padj,2)) %>% 
	mutate(padj = replace_na(padj, 1)) %>% 
	arrange(padj) %>%
	as.data.frame()
rownames(limmaProt) <- limmaProt$UniProt_ID

# Write tables for download
# write.table(limmaRna, "www/mRNA_download.txt", sep="\t", row.names = FALSE)
# write.table(limmaProt, "www/protein_download.txt", sep="\t", row.names = FALSE)

#####
# Misc. hard coded parameters
startGene <- "Q552J0" # random gene here, CHANGE THIS!
pieRadius <- 2.5
loessSpan <- 0.6
tableHeight <- "450px"


#############################
# Define UI
ui <- fluidPage(
	
	tags$head(HTML("<link rel=\"shortcut icon\" type=\"image\" href=\"/dicty_cartoon2.png\"/>")),
	
	title = "Levels of mRNA and protein during Dictyostelium discoideum development",	
	titlePanel(div(HTML("Levels of mRNA and protein during <em>Dictyostelium discoideum</em> development"))),
	
	fluidRow(
		column(3,
					 tags$img(src="dicty_cartoon2.png", alt="logo", width="220px", height="100px")
		),
		
		column(7,
					 p(""),
					 em("N.B."),
					 tags$ul(
					 	tags$li("Points in each scatter plot reflect individual biological replicates"),
					 	tags$li("Significance of each result can be checked using", tags$em("padj")),
					 	tags$li("Imputed values are labelled in gray")
					 ))
		,
		column(2,
					 p(""),
					 downloadButton("report", "Save as pdf")
		)
	),
	
	fluidRow(
		column(5, 
					 plotOutput("mRNAPlot")
		),
		column(5,
					 plotOutput("proteinPlot")
		),
	),
	
	fluidRow(
		tabsetPanel(type = "tabs",
								tabPanel("Select mRNA", value="mrnaTab",
												 tagList("", a("Download table", target="_blank", href="mRNA_download.xlsx")),
												 DTOutput('limmaRnaTab')),
								tabPanel("Select protein", value="proteinTab",
												 tagList("", a("Download table", target="_blank", href="protein_download.xlsx")),
												 DTOutput('limmaProtTab')),
		)
	)
)


#############################  
# Define server
server <- function(input, output, session) {
	
	plotGene <- reactiveVal(startGene)
	
	observeEvent(input$limmaRnaTab_rows_selected, {
		outChar <- as.character(limmaRna[input$limmaRnaTab_rows_selected,"UniProt_ID"])
		if(length(outChar) > 0){
			plotGene(outChar)
		}
	})
	
	observeEvent(input$limmaProtTab_rows_selected, {
		outChar <- as.character(limmaProt[input$limmaProtTab_rows_selected,"UniProt_ID"])
		if(length(outChar) > 0){
			plotGene(outChar)
		}
	})
	
	output$mRNAPlot <- renderPlot({mrnaGgplot()})
	output$proteinPlot <- renderPlot({proteinGgplot()})
	
	mrnaGgplot <- reactive({
		plot.gene(geneId = plotGene(), 
							geneExp = mRNAData[,useRnaSamples],
							geneImpute = mRNADataInputed[,useRnaSamples],
							sampleSheet = sampleSheet[sampleSheet$Sample %in% useRnaSamples,],
							dataAnnot = geneSymbol,
							yaxisUnit = "normalized mRNA expression",
							caption="mRNA",
							span = loessSpan)
	})
	
	proteinGgplot <- reactive({
		plot.gene(geneId = plotGene(), 
							geneExp = proteinData[,useProteinSamples], 
							geneImpute = proteinDataImputed[,useProteinSamples],
							sampleSheet = sampleSheet[sampleSheet$Sample %in% useProteinSamples,],
							dataAnnot = geneSymbol,
							yaxisUnit = "normalized protein expression",
							caption = "Protein",
							span = loessSpan)
	})
	
	output$limmaRnaTab <- DT::renderDataTable(limmaRna, 
																						server = T, 
																						selection = list(mode='single',  selected=which(limmaRna[,1]==startGene)), 
																						rownames = F, 
																						fillContainer = TRUE,
																						options = list(searchHighlight = T, 
																													 fixedHeader=TRUE,
																													 pageLength = nrow(limmaRna),
																													 scrollY = tableHeight,
																													 dom="ft"))
	
	output$limmaProtTab <- DT::renderDataTable(limmaProt, 
																						 server = T, 
																						 selection = list(mode='single',  selected=which(limmaProt[,1]==startGene)), 
																						 rownames = F,
																						 fillContainer = TRUE,
																						 options = list(searchHighlight = T,
																						 							 fixedHeader=TRUE,
																						 							 pageLength = nrow(limmaProt),
																						 							 scrollY = tableHeight,
																						 							 dom="ft"))
	
	mergedPlot <- reactive({ 
		grid.draw(grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=0)),
											 arrangeGrob(mrnaGgplot(), proteinGgplot(),  nrow=2))) 
	})
	
	output$report <- downloadHandler(
		filename = "report.pdf",
		content = function(file) {
			pdf(file=file, width=8, height = 10)
			mergedPlot()
			dev.off()
		}
	)
}

#############################
# Run the application
shinyApp(ui = ui, server = server)