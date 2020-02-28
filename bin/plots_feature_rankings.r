setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)


# barplot of 
barPlotOfRanking <- function(rankingDataFrame, valuesToRank, featuresCol, positionCol, exe_num, x_label, y_label, x_upper_lim) 
{
        library("ggplot2")

        dotSize <- 3
        pdfHeight <- 10 # inches
        pdfWidth <- 20 # inches
        textSize <- 30
        
        mkDirResultsCommand <- "mkdir -p ../results/"
        system(mkDirResultsCommand)
        cat("just run the command: ", mkDirResultsCommand, "\n", sep="")
        
        thisPdfFile <-  paste("../results/", y_label, "_features_", exe_num, ".pdf", sep="")
        pdf(thisPdfFile)
        p <- ggplot(data=rankingDataFrame, aes(x=reorder(featuresCol, -positionCol), y=valuesToRank)) +  geom_bar(stat="identity", fill="steelblue")  + labs(title = paste("Feature ranking on ", y_label, sep=""), y = y_label, x = x_label)        
        
        if ( x_upper_lim !=-1 ) {
            p <- p + coord_flip(ylim = c(0, x_upper_lim))
            cat("\n\n\nciao\n\n\n\n")
        } else {
            p <- p + coord_flip()
        }
        
        
        plot(p)
        cat("saved plot in file ", thisPdfFile, "\n", sep="")
        dev.off()
        
        return(p)
}


fileName <- "../results/LucaOneto_RF_merged_rankings_values.csv"
aggregateRankings <- read.table(file=fileName,head=TRUE,sep=",",stringsAsFactors=FALSE)
print(cat("fileName = ", fileName, "\n", sep=""))

print(cat("colnames(aggregateRankings)\n"))
print(cat(colnames(aggregateRankings)))

num_to_return <- 1
exe_num <- sample(1:upper_num_limit, num_to_return)


FEATURE_RANKING_PLOT_DEPICTION <- TRUE

if (FEATURE_RANKING_PLOT_DEPICTION == TRUE) {
    
        # print(colnames(dd_sorted_IncNodePurity_only))

        mkdirResultsCommand <- "mkdir -p ../results"
        system(mkdirResultsCommand)
        cat("applied command: ", mkdirResultsCommand, "\n", sep="")
          
         upper_num_limit <- 100
         pGini <- barPlotOfRanking(aggregateRankings, aggregateRankings$accuracy_decrease_value, aggregateRankings$accuracy_feature, aggregateRankings$accuracy_pos, exe_num, "features", "accuracy_decrease_value", upper_num_limit)
         
         upper_num_limit <- 100
         accuracyP <- barPlotOfRanking(aggregateRankings, aggregateRankings$Gini_decrease_value, aggregateRankings$Gini_feature, aggregateRankings$Gini_pos, exe_num, "features", "Gini_decrease_value", upper_num_limit)
            
}
        
