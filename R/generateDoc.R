library(rmarkdown)

#rmarkdown::render("chipQC.R",  params = list(file   = '/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChip.csv',name='ALL', consensus=T,summits=250), "github_document")
# /*
#rmarkdown::render("diffBind.R",  "github_document")
#rmarkdown::render("VennPeaks.R", params = list(file   = './test/testUSEQ.xls',
#                                               file2  = './test/testMACSrep1.xls',
#                                               file3  = './test/testMACSrep2.xls',
#                                               outName = 'TEST_VennPeaks.pdf'),
#                                               "github_document")
# */

rmarkdown::render("USEQ.R",  params = list(file='./test/testUSEQ.xls',
                                         outNameBed="TEST_USEQ") 
                                          ,"github_document")
