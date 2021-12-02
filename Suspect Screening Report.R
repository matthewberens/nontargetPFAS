########################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#### Script for reporting suspect screening results ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
########################################################


####### Last edit: 4/1/2019 ############################


# Workflow: First, take the suspect screening .txt file produced with RT of 10 min
# and visualize the data. Then, with SCIEX OS open, determine which compounds
# appear to be hits and create a shortend XIC list with hits and their RTs.
# Alternatively, a shortened XIC list may be provided. In this instance,
# identify likely RTs for the compounds on the shortened XIC list and skip
# the visualization step.
# Run the suspect screening with the shorted XIC list and accurate RTs
# and export that data as a .txt file. Process the second .txt file in 
# this script, starting at row 192, and produce the final excel sheets here.




########################################################
############## Initial Set up ##########################
########################################################

# Set the working directory
setwd("/Users/andrewmaizel/Documents/R working directory/Fluoros/col_exps")

# Load Libraries
library(reshape2) # Provides ability to reshape data
library(ggplot2) # Useful for figures
library(gridExtra) # will allow printing pdfs
library(openxlsx)
library(dplyr)

# Load current XIC list
xic <- read.csv("xic 20180522.csv")


############################################################

# # # # # Produce suspect dashboard for the hits # # # # # 

############################################################



# Load data file 
sus_data <-read.table("20180820 suspect 20190401.txt", 
                               header = T, fill = T, sep = "\t", 
                               stringsAsFactors = F) #Target data export from SCIEX OS



# For visualization, how many of the top X classes (in total peak area) should be displayed?

class_number <-150

# For visualization, how should areas be normalized?

normalize <- "log"  # "max" normalizes to maximum of each analyte
                    # "log" normalizes to a log10 scale,
                    # "no" results in no normalizing

# For visualization, specify the limits on which mass spectral feautures should be included.

ppm_lim <- 10 # Maximum parent mass error (ppm)
iso_lim <- 25 # Maximum  parent isotope error (%) 
area_lim <- 100 # Minimum  parent area
ret_lim <- 5 # Minimum retention time (min)


# Name the .pdf output

plot_name <- "20180820 suspect 20190401.pdf" # What should the final .pdf file be names




############################################################

### Make figures from the suspect screening data file ####

############################################################

# Extract only the mass spectral features which meet the visualization requirements
sustab <- sus_data[which(abs(as.numeric(as.character(sus_data$Mass.Error..ppm.))) < ppm_lim &
                                    as.numeric(as.character(sus_data$Isotope.Ratio.Difference)) < iso_lim &
                                    as.numeric(as.character(sus_data$Area)) > area_lim &
                                    as.numeric(as.character(sus_data$Retention.Time)) > ret_lim
),]

# Remove non-essential columns from the datatable
sustab <- sustab[colnames(sustab) %in% c("Sample.Index", "Sample.Name", "Component.Name",
                                         "Area", "Component.Group.Name", "Precursor.Mass", 
                                         "Retention.Time", "Mass.Error..ppm.", "Library.Score",
                                         "Isotope.Ratio.Difference"
                                         ) 
                 ]

# Melt the measured data so that each analyte measurement gets its own row
sustab <- melt(sustab, id.vars = c("Sample.Index", "Sample.Name", 
                                   "Component.Name", "Component.Group.Name", 
                                   "Precursor.Mass")) 

# Change the column name of "variable" to "measurement"
colnames(sustab)[colnames(sustab) == "variable"] <- "measurement" 

# Replace instances of N/A or infinity with 0
sustab$value[sustab$value == "N/A" | sustab$value == "Infinity"] <- 0 

# Change the value column to numeric
sustab$value <- as.numeric(sustab$value) 

## Determine the top PFAS classes in each sample according to peak area ##
top <- sustab[sustab$measurement == "Area",]
top <- aggregate(value ~ Component.Group.Name, top, sum)
top <- head(top[order(-top$value),],class_number)
names(top)[2] <- "Total Areas"

# Make a data table that is limited to just the compounds of interest
figtab <- sustab[sustab$Component.Group.Name %in% unique(top$Component.Group.Name),]

# normalize for maximum value within each analyte
if (normalize == "max") {
  figtab_areas <- figtab[figtab$measurement == "Area",]
  figtab_areas$value <- figtab_areas$value / aggregate(figtab_areas$value, 
                                                       by = list(figtab_areas$Component.Name),
                                                       max)[match(figtab_areas$Component.Name, 
                                                                  aggregate(figtab_areas$value,
                                                                            by = list(figtab_areas$Component.Name),max)$Group.1  ),]$x 
  figtab <- rbind(figtab[figtab$measurement != "Area",], figtab_areas)
  rm(figtab_areas)
}

# normalize areas on log scale
if (normalize == "log") {
  figtab_areas <- figtab[figtab$measurement == "Area",]
  figtab_areas$value <- log10(figtab_areas$value)
  figtab <- rbind(figtab[figtab$measurement != "Area",], figtab_areas)
  rm(figtab_areas)
}

facet_names <- c('Area' = "Area", 'Retention.Time' = "Ret. Time", 'Mass.Error..ppm.' = "Mass Error", 
                 'Library.Score' = "Lib. Score" ,'Isotope.Ratio.Difference' = "Iso. Ratio Diff.")


pdf(as.character(plot_name), onefile = TRUE)


for (i in 1:length(unique(top$Component.Group.Name))){
  
  figtab_temp <- figtab[figtab$Component.Group.Name == top$Component.Group.Name[i],] 
  
  figtab_temp$Component.Name <- factor(figtab_temp$Component.Name, 
                                       levels = subset(figtab_temp, 
                                                       !duplicated(Component.Name))[order(subset(figtab_temp, 
                                                                                                 !duplicated(Component.Name))$Precursor.Mass),]$Component.Name)
  
  fig <- ggplot(figtab_temp, aes(x = Sample.Index, y = value, fill = Component.Name)) # create a scatter plot
  fig <- fig + geom_line (size = 1, colour = "black") + geom_point(colour = "black", stroke = 2, size = 2, shape = 21) # Format data markers
  fig <- fig + scale_fill_manual(values = rainbow( length( unique(figtab_temp$Component.Name) ) ), name = "Compounds")
  fig <- fig + theme_bw() # removes background
  fig <- fig + coord_cartesian()  # set axes
  fig <- fig + theme(panel.grid = element_blank())  # removes gridlines
  fig <- fig + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"), 
                     plot.title = element_text(face = "bold", size = 20, hjust=.5),
                     legend.title = element_text(face = "bold", size = 14),
                     legend.text = element_text(size = 12)) # Format axes, title
  fig <- fig + theme(panel.border = element_rect(linetype = "solid", colour = "black", size=2))
  fig <- fig + theme(axis.ticks = element_line(linetype= "solid", colour = "black", size = 1.5),
                     axis.ticks.length = unit(.2, "cm")) 
  fig <- fig + labs(x = "Sample Index", y = element_blank(),
                    title = as.character(top$Component.Group.Name[i])) # label the axes
  print(fig + facet_grid(figtab_temp$measurement~., scales = "free_y", labeller = as_labeller(facet_names)) + theme(strip.background = element_blank(), 
                                                                                                                    strip.text = element_text(face = "bold", size = 12)))
  fig
}
dev.off()



# Print a key that associates the sample index values and the sample names, this will
# help identify blanks
write.csv(top, "top.csv")
write.csv(unique(paste(sustab$Sample.Index, sustab$Sample.Name)), "sampleindex.csv")


############################################################

# # # # Make a datatable that can be shared # # # # #  # # # #

############################################################



# Load second SCIEX OS suspect output 
sus_data_short <-read.table("20180820 suspect 20190401.txt", 
                            header = T, 
                            fill = T, 
                            sep = "\t", 
                            stringsAsFactors = F) 

# Load current XIC list
xic <- read.csv("xic 20180522.csv")

# Specify which samples to analyze, either by sample.index or name
#sus_data_short <- sus_data_short[sus_data_short$Sample.Index %in% c(6:13),] # LIMIT BY INDEX
#sustab_short <- sustab_short[grepl("C13|C14|C15|C16|C17|C18", sustab_short$Sample.Name),] # LIMIT BY NAME

# Edit the sample names to remove the analysis date
sus_data_short$Sample.Name <- substring(sus_data_short$Sample.Name, 10,30)

# Define general limits
area_lim <- 100 # Specify the minimum acceptible parent area
ret_lim <- 5 # Specify the minimum retention time

# Define XIC hit limits
xic_ppm_lim <- 5 # Specify the maximum acceptible parent mass error
xic_iso_lim <- 10 # Specify the maximum acceptible parent isotope error  

# Specify Library hit limits
lib_ppm_lim <- 10 # Specify the maximum acceptible parent mass error
lib_iso_lim <- 20 # Specify the maximum acceptible parent isotope error  
lib_lib_lim <- 70 # Specify the minimum libary hit score


############################################

#### Create "Suspect Data (Full)" Table ####

############################################


# Remove columns that are not required for further processing
sustab_short <- sus_data_short[colnames(sus_data_short) %in% 
                               c( "Component.Group.Name", "Component.Name" , "Precursor.Mass", 
                                  "Sample.Name", "Area", "Retention.Time", "Mass.Error..ppm.", 
                                  "Isotope.Ratio.Difference", "Library.Score", "Library.Hit", 
                                  "Formula"
                                  )
                               ]

# Convert value columns to be numeric, with ajust infitite and NA values to 0
sustab_short$Area <- as.numeric(sustab_short$Area)
sustab_short$Retention.Time <- as.numeric(as.character(sustab_short$Retention.Time))
sustab_short$Mass.Error..ppm. <- as.numeric(as.character(sustab_short$Mass.Error..ppm.))
sustab_short$Isotope.Ratio.Difference <- as.numeric(as.character(sustab_short$Isotope.Ratio.Difference))
sustab_short$Library.Score <- as.numeric(as.character(sustab_short$Library.Score))

is.na(sustab_short) <- sapply(sustab_short, is.infinite)

sustab_short$Area[is.na(sustab_short$Area)] <- 0
sustab_short$Retention.Time[is.na(sustab_short$Retention.Time)] <- 0
sustab_short$Mass.Error..ppm.[is.na(sustab_short$Mass.Error..ppm.)] <- 0
sustab_short$Isotope.Ratio.Difference[is.na(sustab_short$Isotope.Ratio.Difference)] <- 0
sustab_short$Library.Score[is.na(sustab_short$Library.Score)] <- 0

# Identify rows as either "Lib", "XIC", or "None"

sustab_short$Hit.Quality <- "None"

sustab_short$Hit.Quality[sustab_short$Area > area_lim &
                           abs(sustab_short$Mass.Error..ppm.) < lib_ppm_lim &
                           sustab_short$Isotope.Ratio.Difference < lib_iso_lim &
                           sustab_short$Library.Score > lib_lib_lim] <- "Lib"

sustab_short$Hit.Quality[sustab_short$Area > area_lim &
                           abs(sustab_short$Mass.Error..ppm.) < xic_ppm_lim &
                           sustab_short$Isotope.Ratio.Difference < xic_iso_lim &
                           sustab_short$Library.Score < lib_lib_lim] <- "XIC"

# Remove rows for compounds that have no library or xic hits

sustab_short <- sustab_short[sustab_short$Component.Name %in% 
                               unique(sustab_short$Component.Name[sustab_short$Hit.Quality %in% c("XIC", "Lib")]),]

sustab_short$In.Lib <- xic$Library[match(sustab_short$Component.Name,xic$Individual.Acronym)]

sustab_short$Super.Class <- xic$Super.Class.2[match(sustab_short$Component.Name,xic$Individual.Acronym)]

sustab_short <- sustab_short[order(sustab_short$Precursor.Mass),]
sustab_short <- sustab_short[order(sustab_short$Component.Group.Name),]
sustab_short <- sustab_short[order(sustab_short$Super.Class),]

# Make tables of average values and standard deviations

sustab_summary <- sustab_short[sustab_short$Hit.Quality != "None",]

sustab_summary <- sustab_summary[,c("Component.Name","Area", "Retention.Time", "Mass.Error..ppm.", "Library.Score","Isotope.Ratio.Difference")]

averages <- aggregate(data = sustab_summary, .~Component.Name,  mean)

stdevs <- aggregate(data = sustab_summary, .~Component.Name,  FUN = sd )

# Determine the most common library match for each compound in the averages table

sustab_lib <- sustab_short[,c("Component.Name","Library.Hit")]
sustab_lib <- sustab_lib[sustab_lib$Library.Hit != "No Match",]
sustab_lib <- sustab_lib[sustab_lib$Library.Hit !="",]
sustab_lib$count <- 1
sustab_lib <- aggregate(data = sustab_lib, .~Component.Name+Library.Hit,  FUN = sum )
sustab_lib <- aggregate(data = sustab_lib, .~Component.Name+Library.Hit, FUN = max)

# Add the most commmon hit to the averages page
averages$Library.Hit <- sustab_lib$Library.Hit[match(averages$Component.Name,sustab_lib$Component.Name)]
averages$Library.Hit[is.na(averages$Library.Hit)] <- "No Match"

# Convert the averages page to something that can be given to clients
averages$Area <- paste(as.character(signif(as.numeric(averages$Area), digits = 2)), " ± " ,
                                    as.character(signif(as.numeric(stdevs$Area), digits = 2)))

averages$Retention.Time <- paste(as.character(signif(as.numeric(averages$Retention.Time), digits = 3)), " ± " ,
                       as.character(signif(as.numeric(stdevs$Retention.Time), digits = 1)))

averages$Mass.Error..ppm. <- paste(as.character(signif(as.numeric(averages$Mass.Error..ppm.), digits = 3)), " ± " ,
                                 as.character(signif(as.numeric(stdevs$Mass.Error..ppm.), digits = 2)))

averages$Isotope.Ratio.Difference <- paste(as.character(signif(as.numeric(averages$Isotope.Ratio.Difference), digits = 2)), " ± " ,
                                 as.character(signif(as.numeric(stdevs$Isotope.Ratio.Difference), digits = 2)))

averages$Library.Score <- paste(as.character(signif(as.numeric(averages$Library.Score), digits = 2)), " ± " ,
                                 as.character(signif(as.numeric(stdevs$Library.Score), digits = 2)))

# Re-order the averages sheet so that it is in the same order as sustab_short
averages <- averages[order(match(averages$Component.Name,unique(sustab_short$Component.Name))),]


# Make a compound info page
cmpd.info <- sustab_short[sustab_short$Hit.Quality != "None",]
cmpd.info <- cmpd.info[,c("Component.Group.Name","Component.Name","Precursor.Mass", "Formula")]
cmpd.info <- unique(cmpd.info)
cmpd.info$Full.Name <- xic$Compound.Name[match(cmpd.info$Component.Name,xic$Individual.Acronym)]
cmpd.info$Group.Name <- xic$Class.Name[match(cmpd.info$Component.Group.Name,xic$Class.Acronym)]
cmpd.info$Super.Class <- xic$Super.Class.2[match(cmpd.info$Component.Name,xic$Individual.Acronym)]
cmpd.info$formula <- xic$Formula[match(cmpd.info$Component.Name,xic$Individual.Acronym)]

# Lastly, add a chunk that identifies structural isomers from the XIC list and adds this to the cmpd.info dataframe

cmpd.info$C <- 0
cmpd.info$H <- 0
cmpd.info$O <- 0
cmpd.info$F <- 0
cmpd.info$S <- 0
cmpd.info$N <- 0
cmpd.info$Cl <- 0

# Determine number of fluorines in each compound
cmpd.info$C <- sub(".*C"  ,"", cmpd.info$Formula) # Remove everything up to "F"
cmpd.info$C <- sub("[A-Z].*$", "", cmpd.info$C) # Remove the first character and everything after it
cmpd.info$C[cmpd.info$C == ""] <- 0

cmpd.info$H <- sub(".*H"  ,"", cmpd.info$Formula) # Remove everything up to "F"
cmpd.info$H <- sub("[A-Z].*$", "", cmpd.info$H) # Remove the first character and everything after it
cmpd.info$H[cmpd.info$H == ""] <- 0

cmpd.info$O <- sub(".*O"  ,"", cmpd.info$Formula) # Remove everything up to "F"
cmpd.info$O <- sub("[A-Z].*$", "", cmpd.info$O) # Remove the first character and everything after it
cmpd.info$O[cmpd.info$O == ""] <- 0

cmpd.info$F <- sub(".*F"  ,"", cmpd.info$Formula) # Remove everything up to "F"
cmpd.info$F <- sub("[A-Z].*$", "", cmpd.info$F) # Remove the first character and everything after it
cmpd.info$F[cmpd.info$F == ""] <- 0

cmpd.info$S <- sub(".*S"  ,"", cmpd.info$Formula) # Remove everything up to "F"
cmpd.info$S <- sub("[A-Z].*$", "", cmpd.info$S) # Remove the first character and everything after it
cmpd.info$S[cmpd.info$S == ""] <- 0

cmpd.info$N <- sub(".*N"  ,"", cmpd.info$Formula) # Remove everything up to "F"
cmpd.info$N <- sub("[A-Z].*$", "", cmpd.info$N) # Remove the first character and everything after it
cmpd.info$N[cmpd.info$N == ""] <- 0

cmpd.info$Cl <- sub(".*Cl"  ,"", cmpd.info$Formula) # Remove everything up to "F"
cmpd.info$Cl <- sub("[A-Z].*$", "", cmpd.info$Cl) # Remove the first character and everything after it
cmpd.info$Cl[cmpd.info$Cl == ""] <- 0

cmpd.info$formula.full<- paste ("C", cmpd.info$C, "H", cmpd.info$H, "O", cmpd.info$O,
                            "F", cmpd.info$F, "S", cmpd.info$S, "N", cmpd.info$N, "Cl", cmpd.info$Cl)

cmpd.info <- cmpd.info[ , -which(names(cmpd.info) %in% c("C","H", "O","N", "F", "S", "Cl"))]


# Extract consistent formulas from XIC list
xic$C <- 0
xic$H <- 0
xic$O <- 0
xic$F <- 0
xic$S <- 0
xic$N <- 0
xic$Cl <- 0

# Determine number of fluorines in each compound
xic$C <- sub(".*C"  ,"", xic$Neutral.Mass.Formula) # Remove everything up to "C"
xic$C <- sub("[A-Z].*$", "", xic$C) # Remove the first character and everything after it
xic$C[xic$C == ""] <- 0

xic$H <- sub(".*H"  ,"", xic$Neutral.Mass.Formula) # Remove everything up to "H"
xic$H <- sub("[A-Z].*$", "", xic$H) # Remove the first character and everything after it
xic$H[xic$H == ""] <- 0

xic$O <- sub(".*O"  ,"", xic$Neutral.Mass.Formula) # Remove everything up to "O"
xic$O <- sub("[A-Z].*$", "", xic$O) # Remove the first character and everything after it
xic$O[xic$O == ""] <- 0

xic$F <- sub(".*F"  ,"", xic$Neutral.Mass.Formula) # Remove everything up to "F"
xic$F <- sub("[A-Z].*$", "", xic$F) # Remove the first character and everything after it
xic$F[xic$F == ""] <- 0

xic$S <- sub(".*S"  ,"", xic$Neutral.Mass.Formula) # Remove everything up to "S"
xic$S <- sub("[A-Z].*$", "", xic$S) # Remove the first character and everything after it
xic$S[xic$S == ""] <- 0

xic$N <- sub(".*N"  ,"", xic$Neutral.Mass.Formula) # Remove everything up to "N"
xic$N <- sub("[A-Z].*$", "", xic$N) # Remove the first character and everything after it
xic$N[xic$N == ""] <- 0

xic$Cl <- sub(".*Cl"  ,"", xic$Neutral.Mass.Formula) # Remove everything up to "Cl"
xic$Cl <- sub("[A-Z].*$", "", xic$Cl) # Remove the first character and everything after it
xic$Cl[xic$Cl == ""] <- 0

xic$formula.full<- paste ("C", xic$C, "H", xic$H, "O", xic$O,
                                "F", xic$F, "S", xic$S, "N", xic$N, "Cl", xic$Cl)

xic <- xic[ , -which(names(xic) %in% c("C","H", "O","N", "F", "S", "Cl"))]

# Now identify compounds with identical formula.full

cmpd.info$isomers <- "None"


for ( i in 1:nrow(cmpd.info)){
xic.temp <- xic[-match(cmpd.info$Component.Name[i], xic$Individual.Acronym),]
cmpd.info$isomers[i] <- paste( xic.temp$Individual.Acronym [which(xic.temp$formula.full %in% cmpd.info$formula.full[i])] , collapse = ", ")
    
}

cmpd.info <- cmpd.info[ , !(names(cmpd.info) %in% c("formula.full"))]

# Re-order the columns to clean up the data for final export

sustab_short <- sustab_short[,c("Super.Class", "Component.Group.Name", "Component.Name", 
                        "Precursor.Mass", "Sample.Name", "Area", "Retention.Time",
                        "Mass.Error..ppm.", "Isotope.Ratio.Difference", "Library.Score",
                        "Library.Hit", "In.Lib", "Hit.Quality"
                        )]

averages$Super.Class <- sustab_short$Super.Class[match(averages$Component.Name, sustab_short$Component.Name)]
averages$Component.Group.Name <- sustab_short$Component.Group.Name[match(averages$Component.Name, sustab_short$Component.Name)]
averages$Precursor.Mass<- sustab_short$Precursor.Mass[match(averages$Component.Name, sustab_short$Component.Name)]

averages <- averages[,c("Super.Class", "Component.Group.Name", "Component.Name", 
                                "Precursor.Mass", "Area", "Retention.Time",
                                "Mass.Error..ppm.", "Isotope.Ratio.Difference", "Library.Score",
                                "Library.Hit"
)]



#####################################################

#### Export an excel sheet with the suspect data ####

#####################################################

write.xlsx(list("Suspect Data" = sustab_short[sustab_short$Hit.Quality != "None",],
                "Suspect Data (Full)" = sustab_short,
                "Averages" = averages, 
                "Compound Info" = cmpd.info), 
           file = "20190914_20190907_neg_JacobsThermal_semiquant.xlsx")


