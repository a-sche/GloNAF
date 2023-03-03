#library("devtools")
#install_github("https://github.com/BioGeoMacro/GIFT")

library("GIFT")

api <- "https://gift.uni-goettingen.de/api/extended/" # ask Patrick for resticted api
GIFT_version <- "beta"

# Bibliographic references in GIFT
ref <- GIFT_references(GIFT_version = GIFT_version, api = api)
head(ref)


### Get GIFT checklist data

# Checklists meta data
checklists_meta <- GIFT_checklist(
  taxon_name = "Spermatophyta",
  complete_taxon = TRUE,
  floristic_group = "native",
  complete_floristic = TRUE,
  geo_type = "All",
  ref_excluded = NULL,
  suit_geo = TRUE,
  shp = NULL, coordinates = NULL, overlap = "centroid_inside",
  remove_overlap = FALSE,
  namesmatched = TRUE,
  list_set_only = TRUE,
  GIFT_version = GIFT_version,
  api = api)

head(ref)


# Identify restricted references we don't want to use
restricted <- ref[which(ref$restricted==1 & ref$ref_ID %in% checklists_meta$lists$ref_ID),]

# For example exclude everything restricted but WCVP
restricted_exclude <- restricted$ref_ID[which(!restricted$ref_ID %in% c(10647))]


# Checklists GIFT
checklists <- GIFT_checklist(
  taxon_name = "Spermatophyta",
  complete_taxon = TRUE,
  floristic_group = "native",
  complete_floristic = TRUE,
  geo_type = "All",
  ref_excluded = restricted_exclude,
  suit_geo = TRUE,
  shp = NULL, coordinates = NULL, overlap = "centroid_inside",
  remove_overlap = FALSE,
  namesmatched = TRUE,
  list_set_only = FALSE,
  GIFT_version = GIFT_version,
  api = api)


# Load GloNAF Checklists
load("GloNAF_WCVP_version20230125/wcvp.accepted.taxa_glonaf.region_20230125.RData")
# Will change to combi_naturalized.RData

# The current version of GIFT only has the OBJSICID of the GloNAF shapefile included 
# I need to load the shapefile
glonaf_shape <- sf::st_read("regions_2020-10-28.shp")
head(glonaf_shape)
names(glonaf_shape)[1] <- "glonaf_ID"

# load GIFT - GloNAF overlap
overlap <- GIFT_overlap(resource = "glonaf", GIFT_version = GIFT_version,
                          api = api)


# Chose regions with sufficient overlap
overlap <- dplyr::left_join(overlap, sf::st_drop_geometry(glonaf_shape[,c(1:2)]))

overlap <- overlap[which(overlap$IDregion %in% colnames(wcvp.accepted.taxa_glonaf.region)),]
overlap <- overlap[which(overlap$overlap12 > 0.98 & overlap$overlap21 > 0.98),]

which(duplicated(overlap$entity_ID))
which(duplicated(overlap$glonaf_ID))

overlap <- dplyr::group_by(overlap, glonaf_ID)
overlap <- dplyr::summarize(overlap, entity_ID = entity_ID[which.max(overlap12 + overlap21)], 
                         overlap12 = overlap12[which.max(overlap12 + overlap21)], 
                         overlap21 = overlap21[which.max(overlap12 + overlap21)], 
                         IDregion = IDregion[which.max(overlap12 + overlap21)])
overlap <- dplyr::ungroup(overlap)

# geodata <- GIFT_shape(checklists_centroid_inside[[1]]$entity_ID, GIFT_version = "beta", api = api)


# subset GIFT checklists

checklists[[1]] <- checklists[[1]][which(checklists[[1]]$entity_ID %in% overlap$entity_ID),]
checklists[[2]] <- checklists[[2]][which(checklists[[2]]$entity_ID %in% overlap$entity_ID),]


# join GloNAF lists

glonaf <- tidyr::pivot_longer(wcvp.accepted.taxa_glonaf.region, cols = 2:922, names_to = "IDregion")
glonaf <- glonaf[which(glonaf$value == 1),]

glonaf <- dplyr::right_join(glonaf, overlap[,c("entity_ID","IDregion")])


names(checklists[[2]])
names(glonaf)

checklists[[2]] <- unique(checklists[[2]][c("entity_ID", "work_species")])
checklists[[2]]$native <- 1

glonaf <- unique(glonaf[c("entity_ID", "wcvp_binomial.accepted")])
names(glonaf) <- c("entity_ID", "work_species")
glonaf$native <- 0

checklists[[2]] <- rbind(checklists[[2]], glonaf)


# Remove overlapping regions
length(unique(checklists[[1]]$entity_ID))
nonoverlapping <- unique(GIFT_no_overlap(checklists[[1]]$entity_ID, area_threshold_island = 0, area_threshold_mainland = 100, overlap_threshold = 0.1))
length(nonoverlapping)

checklists[[1]] <- checklists[[1]][which(checklists[[1]]$entity_ID %in% nonoverlapping),]
checklists[[2]] <- checklists[[2]][which(checklists[[2]]$entity_ID %in% nonoverlapping),]

save(checklists, file = "checklists_GIFT_glonaf.RData")
