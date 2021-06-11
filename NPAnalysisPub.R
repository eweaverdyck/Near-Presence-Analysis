# Title: Near Presence Analysis
# Author Details: Author: Eli Weaverdyck, Contact details: eli.weaverdyck@geschichte.uni-freiburg.de
# Script and data info: This script performs near presence analysis to detect statistically significant clusters of chronological material
# Data consist of survey units with presence/absence data about material from multiple periods
# Data were collected during the Molyvoti, Thrace Archaeology Project in 2015
# Copyright statement: This script is released under the GNU General Public License, version 3.

setwd("C:/Users/eweav/Documents/ArcGIS/Projects/NPAnalysisPub")
library(sf)

# Import data and define periods
chron<-read.csv("ArcGIS_Data/chron.csv",header=TRUE)
units<-st_read("ArcGIS_Data/units.shp", fid_column_name= "FID") # Keep the FID field in the shapefile
near<-read.csv("ArcGIS_Data/near.csv")
near$WT<-1/(near$NEAR_DIST+1) # Add a weight field to the near table of 1/distance+1
periods<- c("Prehist","IA","Arch","Clas","Hell","Cl_He","He_ER","ER","MR","LR","LR_EB","Roma","EB","Byz","Otto")

# Specify number of permutations for comparison distributions and cutoff value
perms<-1000
cut<-0.05

# Create two shapefiles to store results, one for detailed results and one for easily displayed results
np.det<-units
np.disp<-units

# Combine units, chron, and near to get chronological data for all neighbor units
units.chron<-as.data.frame(merge(units, chron, by.x="UnitID",by.y="Survey_Uni")) # Merging units and chron by unit ID associates chronological data with the proper FID
near.units.chron<-merge(near, units.chron, by.x="NEAR_FID",by.y="FID")

# Perform near presence analysis and generate displayable results for each period
for (p in periods){
  
  # Calculate observed near presence. This outputs a data frame with columns "Group.1" (IN_FID) and "x" (NP)
  obs.np<-stats::aggregate((near.units.chron$WT*near.units.chron[,p]), by=list(near.units.chron$IN_FID), FUN=mean)
  
  # Generate a comparison distribution of NP scores with material randomly distributed throughout units
    # Create permutations of the presence/absence data, add an FID field, and merge to the near table as above
    perm.chron<-as.data.frame(replicate(perms,sample(x=as.vector(units.chron[,p]),length(units.chron$FID))))
    perm.chron$FID<-units$FID
    near.perm.chron<-merge(near,perm.chron,by.x="NEAR_FID",by.y="FID")
    
    # Calculate NP scores for each permutation
    perm.np<-stats::aggregate((near.perm.chron$WT*near.perm.chron[,paste("V",1:perms,sep="")]),by=list(near.perm.chron$IN_FID),FUN=mean)
  
  # Determine whether the observed NP is high or low as compared to the distribution of permuted NPs.
    all.np<-merge(obs.np,perm.np,by="Group.1")
    row.names(all.np)<-all.np$Group.1 # Make FIDs into row names so that subsequent functions generate named vectors
    perm.great<-apply(X=all.np[grep("V",names(all.np))]>all.np$x, MARGIN=1, FUN=sum)
    perm.less<-apply(X=all.np[grep("V",names(all.np))]<all.np$x, MARGIN=1, FUN=sum)
    
    # Calculate logical vectors determining whether the observed NP is high or low, based on your cutoff value and number of permutations.
    # Low NP means most permutation NPs are greater than the observed NP, whereas high NP means most are lesser
    Hi.NP<-perm.less>perms-(cut*perms)
    Lo.NP<-perm.great>perms-(cut*perms)
    
    # Merge these vectors using names (i.e., FIDs from units), rename the columns, and merge with chronological data
    HiLo.NP<-merge(Hi.NP,Lo.NP,by=0)
    colnames(HiLo.NP)<-c("FID","HiNP","LoNP")
    chron.HiLo.NP<-merge(units.chron[c("FID",p)],HiLo.NP,by="FID")
    
    # Create text field with interpretable, displayable results
    chron.HiLo.NP[,paste(p,"_NP",sep="")]<-ifelse(
      chron.HiLo.NP[,p]==1 & chron.HiLo.NP$HiNP==TRUE, 
        "Present, High NP",
        ifelse(chron.HiLo.NP[,p]==1 & chron.HiLo.NP$HiNP==FALSE & chron.HiLo.NP$LoNP==FALSE, 
               "Present, Moderate NP",
               ifelse(chron.HiLo.NP[,p]==1 & chron.HiLo.NP$LoNP==TRUE,
                      "Present, Low NP",
                      ifelse(chron.HiLo.NP[,p]==0 & chron.HiLo.NP$HiNP==TRUE,
                             "Absent, High NP",
                             ifelse(chron.HiLo.NP[,p]==0 & chron.HiLo.NP$HiNP==FALSE & chron.HiLo.NP$LoNP==FALSE,
                                    "Absent, Moderate NP",
                                    ifelse(chron.HiLo.NP[,p]==0 & chron.HiLo.NP$LoNP==TRUE,
                                           "Absent, Low NP","ERROR")))))
    )
  
   
    
    # Merge the results into the units data sets created to store them.
    np.disp<-merge(np.disp, chron.HiLo.NP[, c("FID", paste(p, "_NP", sep=""))], by.x="FID", by.y="FID")
    row.names(np.disp)<-np.disp$FID
    
    # The detailed results require some renaming
    colnames(obs.np)<-c("FID",paste(p,"_np",sep=""))
    ranks<-merge(perm.great,perm.less,by=0)
    colnames(ranks)<-c("FID",paste(p,"_prmgrtr",sep=""),paste(p,"_prmlssr",sep=""))
    np.det<-merge(np.det,obs.np,by="FID")
    np.det<-merge(np.det,ranks,by="FID")
}

# Determine which units are either cluster cores (Present, High NP) or contribute to the high NP of the cores
np.nogeo<-st_drop_geometry(np.disp) # makes subsetting easier
for (p in periods){
  cores<-np.nogeo[np.nogeo[, paste(p, "_NP", sep="")]=="Present, High NP", "FID"] # Vector of FIDs of units with present, high NP
  for (f in np.disp$FID){
    np.disp[f, paste(p, "_mem", sep="")]<-np.nogeo[f, paste(p, "_NP", sep="")]=="Present, High NP" |
      (sum(near.units.chron[near.units.chron$NEAR_FID == f,p])>0 &
      any(near.units.chron[near.units.chron$NEAR_FID==f,"IN_FID"] %in% cores))
  }
}

# Write the results to shapefiles
st_write(np.disp,"npdisp.shp")
st_write(np.det, "npdet.shp")
