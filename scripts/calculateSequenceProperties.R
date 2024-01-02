suppressWarnings(library(Peptides))
suppressWarnings(library(seqinr))
suppressWarnings(library(stringr))
options(warn=2)
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]

db<-data.frame()
fa<-read.fasta(input,as.string=TRUE)


n=0
for (i in fa){
#	print(getName(i))
	n=n+1
	db[n,'Name']=getName(i)
	db[n,'Length']=lengthpep(seq = i[1])
	db[n,'Number of cytosines']=str_count(i[1],"c")
	db[n,'Number of acidic AA']=str_count(i[1],"d")+str_count(i[1],"e")
	db[n,'Number of basic AA']=str_count(i[1],"k")+str_count(i[1],"r")+str_count(i[1],"h")
	
	db[n,'Number of aromatic AA']=str_count(i[1],"f")+str_count(i[1],"y")+str_count(i[1],"w")
	db[n,'Number of polar and uncharged AA']=str_count(i[1],"s")+str_count(i[1],"t")+str_count(i[1],"c")+str_count(i[1],"p")+str_count(i[1],"n")+str_count(i[1],"q")
	db[n,'Number of nonpolar and aliphatic AA']=str_count(i[1],"g")+str_count(i[1],"a")+str_count(i[1],"v")+str_count(i[1],"l")+str_count(i[1],"m")+str_count(i[1],"i")
	db[n,'Molecular Weight']=mw(seq=i[1])
	db[n,'Net Charge']=charge(seq=i[1], pH = 7, pKscale = "EMBOSS")
	db[n,'Isoelectric Point (pI)']=pI(seq=i[1], pKscale = "EMBOSS")
	db[n,'Aliphatic Index']=pI(seq=i[1])
	db[n,'Instability Index']=instaIndex(seq = i[1])
	db[n,'Boman Index']=boman(seq = i[1])
	db[n,'Hydrophobicity']=hydrophobicity(seq = i[1], scale = "Eisenberg")
	db[n,'Hydrophobic Moment']=hmoment(seq = i[1], angle = 100, window = 11)
}	
write.csv(db,output)
