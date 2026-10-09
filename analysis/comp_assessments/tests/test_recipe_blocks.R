source('analysis/comp_assessments/R/run_assessment.R')
db <- read_database()
translate <- function(id) {
 e <- new.env(parent = environment())
 sys.source(file.path('analysis/comp_assessments/scripts/translation/stocks',paste0(id,'.R')),e)
 e$translate_stock(read_assessment(id,db))$obs
}
expect_block <- function(actual, expected) stopifnot(identical(as.character(actual),as.character(expected)))

obs <- translate('afsc_pollock_goa_2024')
expect_block(obs$index$q_age_block,paste0(2*((obs$index$age-1)%/%2)+1,'-',2*((obs$index$age-1)%/%2)+2))
adfg <- 'ADF&G crab/groundfish trawl'
expect_block(obs$index$adfg_year,ifelse(obs$index$survey==adfg,obs$index$year,min(obs$index$year[obs$index$survey==adfg])))
obs <- translate('afsc_pollock_ebs_2024')
expect_block(obs$index$q_age_block,ifelse(obs$index$age>=9,'9+',obs$index$age))
obs <- translate('dfo_cod_2j3kl_2025')
rv <- obs$index$survey == 'DFO fall RV survey'
expected <- paste('age',paste0(2*((obs$index$age-2)%/%2)+2,'-',pmin(2*((obs$index$age-2)%/%2)+3,14)))
expected[obs$index$age>=12] <- 'age 12-14'
expected[!rv] <- 'age 2-3'
expect_block(obs$index$smith_sound_q_key,expected)
obs <- translate('ices_bluewhiting_northeast_atlantic_2026')
expect_block(obs$catch$sd_group,ifelse(obs$catch$age<=2,obs$catch$age,ifelse(obs$catch$age<=8,'3-8','9-10')))
expect_block(obs$index$sd_group,ifelse(obs$index$age<=3,obs$index$age,ifelse(obs$index$age<=6,'4-6','7-8')))
expect_block(obs$index$q_key,paste(obs$index$survey,ifelse(obs$index$age<=4,obs$index$age,'5-8'),sep='.'))
obs <- translate('ices_herring_north_sea_2026')
expect_block(obs$catch$sd_block,ifelse(obs$catch$age<=1,'ages_0_1',ifelse(obs$catch$age<=6,'ages_2_6','ages_7_8')))
s <- obs$index$survey; a <- obs$index$age
expect_block(obs$index$sd_block,ifelse(s=='HERAS',paste0('HERAS_',ifelse(a<=3,a,ifelse(a<=6,'4_6','7_8'))),ifelse(s=='IBTS-Q3',paste0('IBTS-Q3_',ifelse(a<=1,'0_1','2_5')),s)))
expect_block(obs$index$q_key,ifelse(s=='HERAS',ifelse(a<=2,'HERAS_1_2','HERAS_3_8'),ifelse(s=='IBTS-Q3',paste0(s,'_',a),s)))
obs <- translate('ices_herring_western_baltic_2026')
expect_block(obs$catch$sd_block,ifelse(obs$catch$age==0,'sam_sd_5',ifelse(obs$catch$age==1,'sam_sd_6','sam_sd_0')))
s <- obs$index$survey; a <- obs$index$age
expect_block(obs$index$q_key,ifelse(s=='HERAS',ifelse(a==2,'sam_q_0','sam_q_1'),ifelse(s=='GERAS',ifelse(a<=2,'sam_q_2','sam_q_7'),ifelse(s=='N20','sam_q_3',paste0('sam_q_',a+1)))))
obs <- translate('ices_norway_pout_north_sea_2026_benchmark')
expect_block(obs$catch$sd_block,ifelse(obs$catch$age==0,'age0',ifelse(obs$catch$age<=2,'age1_2','age3plus')))
obs <- translate('ices_saithe_north_sea_2026')
expect_block(obs$catch$sd_block,ifelse(obs$catch$age==3,'age3',ifelse(obs$catch$age<=5,'age4_5','age6_plus')))
obs <- translate('ices_sprat_baltic_2026')
expect_block(obs$index$q_key,paste(obs$index$survey,ifelse(obs$index$age>=6,'6-8',obs$index$age),sep='.'))
cat('Recipe age and year block boundaries are unchanged.\n')

