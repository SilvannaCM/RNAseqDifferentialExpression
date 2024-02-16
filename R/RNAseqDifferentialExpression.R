#hola
# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)

# Which data do you want to use? Let's use the sailfish counts.
# browseURL("http://dx.doi.org/10.6084/m9.figshare.1601975")
# countDataURL = "http://files.figshare.com/2439061/GSE37704_featurecounts.csv"
countDataURL = "/Users/silvanacristo/Downloads/GSE37704_featurecounts.csv"

# Import countdata
countData = read.csv(countDataURL, row.names=1) %>%
  dplyr::select(-length) %>%
  as.matrix()

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)

##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000198888     17528     23007     30241     24418     29152
## ENSG00000198763     21264     26720     35550     28878     32416
## ENSG00000198804    130975    151207    195514    178130    196727
## ENSG00000198712     49769     61906     78608     66478     69758
## ENSG00000228253      9304     11160     12830     12608     13041
## ENSG00000198899     45401     51260     66851     63433     66123
##                 SRR493371
## ENSG00000198888     34416
## ENSG00000198763     38422
## ENSG00000198804    244670
## ENSG00000198712     86808
## ENSG00000228253     16063
## ENSG00000198899     79215

# Import metadata
colData = read.csv("/Users/silvanacristo/Downloads/GSE37704_metadata.csv", row.names=1)
colData

##               condition
## SRR493366 control_sirna
## SRR493367 control_sirna
## SRR493368 control_sirna
## SRR493369      hoxa1_kd
## SRR493370      hoxa1_kd
## SRR493371      hoxa1_kd

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
dds = DESeq(dds)
dds

## class: DESeqDataSet
## dim: 16755 6
## metadata(0):
## assays(3): counts mu cooks
## rownames(16755): ENSG00000198888 ENSG00000198763 ...
##   ENSG00000267795 ENSG00000165795
## rowRanges metadata column names(27): baseMean baseVar ... deviance
##   maxCooks
## colnames(6): SRR493366 SRR493367 ... SRR493370 SRR493371
## colData names(2): condition sizeFactor

res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
res = res[order(res$pvalue),]
summary(res)

##
## out of 16755 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4193, 25%
## LFC < 0 (down)   : 4286, 26%
## outliers [1]     : 22, 0.13%
## low counts [2]   : 1299, 7.8%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results

# Instalar paquetes antes de usarlos, sin embargo no estan disponibles asi que usare una funcion extra.

# Verificar si BiocManager está instalado
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  print("BiocManager no está instalado. Instalando...")
  install.packages("BiocManager")
} else {
  print("BiocManager está instalado.")
}

# Funcion extra
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

# Verificar si los paquetes están instalados
installed_packages <- installed.packages()
if ("AnnotationDbi" %in% rownames(installed_packages)) {
  print("El paquete AnnotationDbi está instalado.")
} else {
  print("El paquete AnnotationDbi no está instalado.")
}

if ("org.Hs.eg.db" %in% rownames(installed_packages)) {
  print("El paquete org.Hs.eg.db está instalado.")
} else {
  print("El paquete org.Hs.eg.db no está instalado.")
}

#Sigue el tutorial
columns(org.Hs.eg.db)

##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"
##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"
## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"
## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"
## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
## [25] "UNIGENE"      "UNIPROT"

res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

head(res, 10)

## log2 fold change (MAP): condition hoxa1_kd vs control_sirna
## Wald test p-value: condition hoxa1_kd vs control_sirna
## DataFrame with 10 rows and 9 columns
##                  baseMean log2FoldChange      lfcSE      stat    pvalue
##
## ENSG00000148773  1885.344      -3.172502 0.07868572 -40.31865         0
## ENSG00000138623  2939.936      -2.418238 0.05889229 -41.06205         0
## ENSG00000104368 13601.963       2.016802 0.05249643  38.41789         0
## ENSG00000124766  2692.200       2.379545 0.06193654  38.41908         0
## ENSG00000122861 35889.413       2.224779 0.05258658  42.30697         0
## ENSG00000116016  4558.157      -1.885339 0.04258766 -44.26961         0
## ENSG00000164251  2404.103       3.325196 0.07021236  47.35912         0
## ENSG00000125257  6187.386       1.943762 0.04259189  45.63692         0
## ENSG00000104321  9334.555       3.186856 0.06227530  51.17367         0
## ENSG00000183508  2110.345       3.190612 0.07488305  42.60794         0
##                      padj      symbol      entrez
##
## ENSG00000148773         0       MKI67        4288
## ENSG00000138623         0      SEMA7A        8482
## ENSG00000104368         0        PLAT        5327
## ENSG00000124766         0        SOX4        6659
## ENSG00000122861         0        PLAU        5328
## ENSG00000116016         0       EPAS1        2034
## ENSG00000164251         0       F2RL1        2150
## ENSG00000125257         0       ABCC4       10257
## ENSG00000104321         0       TRPA1        8989
## ENSG00000183508         0      FAM46C       54855
##                                                                               name
##
## ENSG00000148773                                      marker of proliferation Ki-67
## ENSG00000138623 semaphorin 7A, GPI membrane anchor (John Milton Hagen blood group)
## ENSG00000104368                                      plasminogen activator, tissue
## ENSG00000124766                               SRY (sex determining region Y)-box 4
## ENSG00000122861                                   plasminogen activator, urokinase
## ENSG00000116016                                   endothelial PAS domain protein 1
## ENSG00000164251                   coagulation factor II (thrombin) receptor-like 1
## ENSG00000125257            ATP-binding cassette, sub-family C (CFTR/MRP), member 4
## ENSG00000104321 transient receptor potential cation channel, subfamily A, member 1
## ENSG00000183508                       family with sequence similarity 46, member C

#Instalar pathview
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")

# Verificar si el paquete está instalado
installed_packages <- installed.packages()
if ("pathview" %in% rownames(installed_packages)) {
  print("El paquete pathview está instalado.")
} else {
  print("El paquete pathview no está instalado.")
}

#Instalar gage y gageData
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gage")
BiocManager::install("gageData")

# Verificar si los paquetes están instalados
installed_packages <- installed.packages()
if ("gage" %in% rownames(installed_packages)) {
  print("El paquete gage está instalado.")
} else {
  print("El paquete gage no está instalado.")
}

if ("gageData" %in% rownames(installed_packages)) {
  print("El paquete gageData está instalado.")
} else {
  print("El paquete gageData no está instalado.")
}

#Sigue el tutorial
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

## $`hsa00232 Caffeine metabolism`
## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"
##
## $`hsa00983 Drug metabolism - other enzymes`
##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"
##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"
## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733"
## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600"
## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"
## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"
## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"
## [50] "8833"   "9"      "978"
##
## $`hsa00230 Purine metabolism`
##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"
##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164"
##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649"
##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115"
##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"
##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"
##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"
##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"
##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"
##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"
##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"
##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"
##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"
## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"
## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"
## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703"
## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985"
## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618"
## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034"
## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"
## [162] "9615"

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

##      4288      8482      5327      6659      5328      2034
## -3.172502 -2.418238  2.016802  2.379545  2.224779 -1.885339

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

## $greater
##                                          p.geomean stat.mean        p.val
## hsa04142 Lysosome                     0.0002630657  3.517890 0.0002630657
## hsa04640 Hematopoietic cell lineage   0.0017919390  2.976432 0.0017919390
## hsa04630 Jak-STAT signaling pathway   0.0048980977  2.604390 0.0048980977
## hsa00140 Steroid hormone biosynthesis 0.0051115493  2.636206 0.0051115493
## hsa04062 Chemokine signaling pathway  0.0125582961  2.250765 0.0125582961
## hsa00511 Other glycan degradation     0.0223819919  2.104311 0.0223819919
##                                            q.val set.size         exp1
## hsa04142 Lysosome                     0.04261664      116 0.0002630657
## hsa04640 Hematopoietic cell lineage   0.14514706       61 0.0017919390
## hsa04630 Jak-STAT signaling pathway   0.20701775      119 0.0048980977
## hsa00140 Steroid hormone biosynthesis 0.20701775       39 0.0051115493
## hsa04062 Chemokine signaling pathway  0.40688879      156 0.0125582961
## hsa00511 Other glycan degradation     0.49956506       15 0.0223819919
##
## $less
##                                      p.geomean stat.mean        p.val
## hsa04110 Cell cycle               2.165725e-06 -4.722301 2.165725e-06
## hsa03030 DNA replication          3.807440e-06 -4.835336 3.807440e-06
## hsa04114 Oocyte meiosis           1.109869e-04 -3.767561 1.109869e-04
## hsa03013 RNA transport            1.181787e-03 -3.071947 1.181787e-03
## hsa03440 Homologous recombination 1.197124e-03 -3.190747 1.197124e-03
## hsa00240 Pyrimidine metabolism    1.570318e-03 -2.992059 1.570318e-03
##                                          q.val set.size         exp1
## hsa04110 Cell cycle               0.0003084027      121 2.165725e-06
## hsa03030 DNA replication          0.0003084027       36 3.807440e-06
## hsa04114 Oocyte meiosis           0.0059932916      101 1.109869e-04
## hsa03013 RNA transport            0.0387868193      145 1.181787e-03
## hsa03440 Homologous recombination 0.0387868193       28 1.197124e-03
## hsa00240 Pyrimidine metabolism    0.0423985796       96 1.570318e-03
##
## $stats
##                                       stat.mean     exp1
## hsa04142 Lysosome                      3.517890 3.517890
## hsa04640 Hematopoietic cell lineage    2.976432 2.976432
## hsa04630 Jak-STAT signaling pathway    2.604390 2.604390
## hsa00140 Steroid hormone biosynthesis  2.636206 2.636206
## hsa04062 Chemokine signaling pathway   2.250765 2.250765
## hsa00511 Other glycan degradation      2.104311 2.104311

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=5) %>%
  .$id %>%
  as.character()
keggrespathways

## [1] "hsa04142 Lysosome"
## [2] "hsa04640 Hematopoietic cell lineage"
## [3] "hsa04630 Jak-STAT signaling pathway"
## [4] "hsa00140 Steroid hormone biosynthesis"
## [5] "hsa04062 Chemokine signaling pathway"

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

## [1] "hsa04142" "hsa04640" "hsa04630" "hsa00140" "hsa04062"

#Crear los plots
# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

