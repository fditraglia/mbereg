# Oreopoulos (2006, AER) ===========================================
download.file('https://www.aeaweb.org/aer/data/mar06_data_20030407.zip',
              './data-raw/Oreopoulos.zip')
unzip('./data-raw/Oreopoulos.zip',
      files = 'uk/combined general household survey.dta',
      exdir = './data-raw')
system('mv ./data-raw/uk/combined\\ general\\ household\\ survey.dta ./data-raw')
system('rm -r ./data-raw/uk')
system('rm ./data-raw/Oreopoulos.zip')
oreopoulos <- haven::read_dta('./data-raw/combined general household survey.dta')

# The following is an *exact* translation of the STATA code in
#     "make cell means data 2008 08 21.do"
# into R. This do-file generates additional variables used in the analysis
# and subsets the data. I do not include the step in which Oreopoulus computes
# cell means, since we will not use these in our analysis.

oreopoulos$yearat14 <- with(oreopoulos, yobirth + 14)
oreopoulos$ageat47 <- with(oreopoulos, 47 - yobirth)
oreopoulos$ageat57 <- with(oreopoulos, 57 - yobirth)
oreopoulos$ageat72 <-  with(oreopoulos, 72 - yobirth)

oreopoulos$drop14 <- with(oreopoulos, ifelse(nireland == 0, ageat47 >= 15,
                                             ageat57 >= 15))

oreopoulos$lhinc <- with(oreopoulos, log(hinc))
oreopoulos$learn <- with(oreopoulos, log(earn))
oreopoulos$linc <- with(oreopoulos, log(inc))
oreopoulos$drop15 <- with(oreopoulos, drop14 == 0)
oreopoulos$drop16 <- with(oreopoulos, ageat72 < 15)

oreopoulos$educb14 <- with(oreopoulos, agelfted <= 14)
oreopoulos$educb15 <- with(oreopoulos, agelfted <= 15)

oreopoulos$badhealth <- with(oreopoulos, (genhealth==3) * (genhealth >= 1 & genhealth <= 3))
oreopoulos$goodhealth <- with(oreopoulos, (genhealth==1) * (genhealth >= 1 & genhealth <= 3))

oreopoulos$dropage <- NA
drop14 <- with(oreopoulos, ((ageat47 >= 15) & (nireland == 0)) |
         ((ageat57 >=15) & (nireland == 1)))
drop15 <- with(oreopoulos, (ageat72 >= 15) & (dropage != 14))
oreopoulos$dropage[drop14] <- 14
oreopoulos$dropage[drop15] <- 15
oreopoulos$dropage[is.na(oreopoulos$dropage)] <- 16
rm(drop14, drop15)

oreopoulos$yearat14_2 <- with(oreopoulos, yearat14^2)
oreopoulos$yearat14_3 <- with(oreopoulos, yearat14^3)
oreopoulos$yearat14_4 <- with(oreopoulos, yearat14^4)
oreopoulos$age2 <- with(oreopoulos, age^2)
oreopoulos$age3 <- with(oreopoulos, age^3)
oreopoulos$age4 <- with(oreopoulos, age^4)

oreopoulos <- subset(oreopoulos, (yearat14 >= 35) & (yearat14 <= 65) & (age <= 64))
oreopoulos <- subset(oreopoulos, (agelfted >= 10) & !is.na(agelfted))
oreopoulos$missing_earn <- with(oreopoulos, is.na(learn))
