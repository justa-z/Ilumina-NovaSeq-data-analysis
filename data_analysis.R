# Wczytanie niezbędnych pakietów 

library(dada2) 
library(ggplot2) 
library(phyloseq) 

# Ustawienie katalogu roboczego 
setwd("D:/. . ./R") 

# Definiowanie ścieżek do plików FASTQ (sekwencje przedniej i tylnej nici) 
fnFs <- list.files(path = "D:/. . ./R", pattern = "_R1_001.fastq", full.names = TRUE) 
fnRs <- list.files(path = "D:/. . ./R", pattern = "_R2_001.fastq", full.names = TRUE) 

# Wizualizacja jakości pierwszego pliku dla obu nici 
plotQualityProfile(fnFs[1]) 
plotQualityProfile(fnRs[1]) 

# Definiowanie ścieżki dla przefiltrowanych plików 
filtered_path <- file.path("D:/. . ./. . .", "filtered" ) 
filtFs <- file.path(filtered_path, paste0(basename(fnFs), "_filtered.fastq")) 
filtRs <- file.path(filtered_path, paste0(basename(fnRs), "_filtered.fastq")) 

# Filtrowanie i przycinanie sekwencji
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240,230),  
                     maxN = 0, maxEE = c(2,6), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread=TRUE) 

# Uczenie modelu błędów dla obu nici 
errF <- learnErrors(filtFs, multithread = TRUE) 
errR <- learnErrors(filtRs, multithread = TRUE) 

# Wizualizacja modelu błędów 
plotErrors(errF, nominalQ = TRUE)

# Dopasowanie sekwencji do modelu błędów 
dadaFs <- dada(filtFs, err = errF, multithread = TRUE) 
dadaRs <- dada(filtRs, err = errR, multithread = TRUE) 

# Łączenie par sekwencji 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE) 
summary(mergers) 

# Tworzenie tabeli ASV (Amplicon Sequence Variant) 
seqtab <- makeSequenceTable(mergers) 
cat("Unikalne sekwencje: ", ncol(seqtab), "\n") 

# Usuwanie chimer 
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "per-sample", multithread = FALSE, verbose = TRUE) 
cat("Sekwencje po filtracji: ", ncol(seqtab.nochim), "\n") 

# Przypisanie taksonomii przy użyciu bazy SILVA 
silva_db <-"D:/. . ./. . ./silva_nr99_v138.1_train_set.fa.gz" 
taxa <- assignTaxonomy(seqtab.nochim, silva_db, multithread = TRUE) 
print(head(taxa)) 
