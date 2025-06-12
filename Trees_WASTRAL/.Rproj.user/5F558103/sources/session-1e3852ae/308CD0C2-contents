# --- Code written by Cassia Bitencourt but also include function provided by Sidonie Bellot ---
# This script is implemented for plotting the summary trees resulting from wASTRAL, and it is written into five steps. [1] First we created vectors which will translate into a block of colours later. 
# This means that with these blocks we will be able to provide different colours for the terminals, in this case, I am splitting them by grades (rauvolfioids and apocynoids), subfamilies (Periplocoideae, 
# Secamonoideae and Asclepiadoideae). Plus the outgroups and rogue taxa. [2] In this step, we provide functions to define the colours per vector, and another function developed by Sidonie Bellot to calculate 
# the local posterior probability and a final function to set up the plot of the phylogenetic tree. The end of this script is for loading the data. The tree that came from the summary method and a CSV file including two 
# columns the current names of the terminals and the new names you would like to use. Here we include subfamily, tribe, subtribe, genus, species and project codification or accession number. Finally, we plot 
# the tree and save it as a PDF.

# --- Load libraries ---
library(ggtree)
library(ape)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(paletteer)  # For custom palettes

# --- Define taxon groups ---
# Define rogue taxa
rogue_taxa <- c("Secamonoideae_notribe_nosubtribe_Goniostemma_acuminatum_PAFTOL025379",
                "Asclepiadoideae_Asclepiadeae_Gonolobinae_Macroscepis_multiflora_PAFTOL031788",
                "Asclepiadoideae_Asclepiadeae_Cynanchinae_Cynanchum_laeve_PAFTOL031815") 

# Define outgroups
outgroups_taxa <- c("Gelsemiaceae_Gelsemieae_Gelsemium_elegans_SRR12009649",
                    "Gentianaceae_Chironieae_Chironia_baccifera_PAFTOL005228",
                    "Gentianaceae_Exaceae_Exacum_exiguum_PAFTOL005229",
                    "Gentianaceae_Helieae_Macrocarpaea_stenophylla_PAFTOL005232",
                    "Gentianaceae_Potalieae_Utania_racemosa_PAFTOL005230",
                    "Gentianaceae_Voyrieae_Voyria_aurantiaca_PAFTOL005233",
                    "Loganiaceae_Antonieae_Antonia_ovata_PAFTOL005873",
                    "Loganiaceae_Loganieae_Geniostoma_borbonicum_PAFTOL006558",
                    "Loganiaceae_Spigelieae_Spigelia_pulchella_PAFTOL005874",
                    "Loganiaceae_Strychneae_Neuburgia_corynocarpa_PAFTOL007837",
                    "Loganiaceae_Strychneae_Strychnos_splendens_PAFTOL005235",
                    "Rubiaceae_Colletoecemateae_Colletoecema_tortistilum_PAFTOL005241",
                    "Rubiaceae_Luculieae_Luculia_pinceana_PAFTOL005309",
                    "Rubiaceae_Ophiorrhizeae_Ophiorrhiza_winkleri_PAFTOL005263",
                    "Rubiaceae_Perameae_Perama_dichotoma_PAFTOL005268",
                    "Rubiaceae_Psychotrieae_Psychotria_borbonica_PAFTOL021587",
                    "Rubiaceae_Putorieae_Plocama_calabrica_PAFTOL021615",
                    "Rubiaceae_Rubieae_Didymaea_alsinoides_PAFTOL032765",
                    "Rubiaceae_Schizocoleeae_Schizocolea_linderi_PAFTOL005887")

# Define rauvolfioids
rauvolfioids_taxa <- c("rauvolfioids_Alstonieae_nosubtribe_Alstonia_macrophylla_PAFTOL004113",
                       "rauvolfioids_Alstonieae_nosubtribe_Alstonia_scholaris_PAFTOL006553",
                       "rauvolfioids_Alstonieae_nosubtribe_Alstonia_spectabilis_PAFTOL007872",
                       "rauvolfioids_Alyxieae_Alyxiinae_Alyxia_buxifolia_PAFTOL004100",
                       "rauvolfioids_Alyxieae_Alyxiinae_Lepinia_solomonensis_PAFTOL006601",
                       "rauvolfioids_Alyxieae_Alyxiinae_Lepinia_taitensis_PAFTOL006555",
                       "rauvolfioids_Alyxieae_Alyxiinae_Pteralyxia_laurifolia_PAFTOL031741",
                       "rauvolfioids_Alyxieae_Condylocarpinae_Chilocarpus_costatus_PAFTOL031728",
                       "rauvolfioids_Alyxieae_Condylocarpinae_Condylocarpon_isthmicum_PAFTOL031721",
                       "rauvolfioids_Alyxieae_Condylocarpinae_Plectaneia_thouarsii_PAFTOL031781",
                       "rauvolfioids_Amsonieae_nosubtribe_Amsonia_hubrichtii_PAFTOL005220",
                       "rauvolfioids_Aspidospermateae_nosubtribe_Aspidosperma_cylindrocarpon_PAFTOL004040",
                       "rauvolfioids_Aspidospermateae_nosubtribe_eissospermum_argenteum_PAFTOL025219",
                       "rauvolfioids_Aspidospermateae_nosubtribe_Haplophyton_cimicidum_PAFTOL031732",
                       "rauvolfioids_Aspidospermateae_nosubtribe_Microplumeria_anomala_PAFTOL031736",
                       "rauvolfioids_Aspidospermateae_nosubtribe_Strempeliopsis_strempelioides_PAFTOL025289",
                       "rauvolfioids_Aspidospermateae_nosubtribe_Vallesia_glabra_PAFTOL031744",
                       "rauvolfioids_Carisseae_nosubtribe_Acokanthera_oblongifolia_PAFTOL025441",
                       "rauvolfioids_Carisseae_nosubtribe_Acokanthera_schimperi_PAFTOL031725",
                       "rauvolfioids_Carisseae_nosubtribe_Carissa_macrocarpa_PAFTOL005221",
                       "rauvolfioids_Carisseae_nosubtribe_Carissa_spinarum_PAFTOL025443",
                       "rauvolfioids_Hunterieae_nosubtribe_Gonioma_kamassi_PAFTOL004050",
                       "rauvolfioids_Hunterieae_nosubtribe_Picralima_nitida_PAFTOL025249",
                       "rauvolfioids_Hunterieae_nosubtribe_Pleiocarpa_rostrata_PAFTOL025251",
                       "rauvolfioids_Melodineae_nosubtribe_Melodinus_forbesii_PAFTOL005890",
                       "rauvolfioids_Melodineae_nosubtribe_Stephanostegia_hildebrandtii_PAFTOL025363",
                       "rauvolfioids_Plumerieae_Plumeriinae_Himatanthus_attenuatus_PAFTOL031692",
                       "rauvolfioids_Plumerieae_Plumeriinae_Mortoniella_pittieri_PAFTOL031794",
                       "rauvolfioids_Plumerieae_Plumeriinae_Plumeria_obtusa_PAFTOL005225",
                       "rauvolfioids_Plumerieae_Thevetiinae_Anechites_nerium_PAFTOL031747",
                       "rauvolfioids_Plumerieae_Thevetiinae_Cerbera_dumicola_GAP022557",
                       "rauvolfioids_Plumerieae_Thevetiinae_Cerberiopsis_obtusifolia_PAFTOL025407",
                       "rauvolfioids_Plumerieae_Thevetiinae_Thevetia_bicornuta_PAFTOL031724",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Ambelania_acida_PAFTOL025325",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Macoubea_guianensis_PAFTOL031734",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Molongum_laxum_PAFTOL025221",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Mucoa_duckei_PAFTOL031810",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Neocouma_ternstroemiacea_PAFTOL031811",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Rhigospira_quadrangularis_PAFTOL025223",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Spongiosperma_grandiflorum_PAFTOL031718",
                       "rauvolfioids_Tabernaemontaneae_Ambelaniinae_Spongiosperma_riparium_PAFTOL025497",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Callichilia_subsessilis_PAFTOL025317",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Calocrater_preussii_PAFTOL025305",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Carvalhoa_campanulata_PAFTOL031753",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Crioceras_dipladeniiflorus_PAFTOL025279",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Ephippiocarpa_orientalis_PAFTOL031764",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Schizozygia_coffeoides_PAFTOL025453",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Tabernaemontana_divaricata_PAFTOL005308",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Tabernanthe_iboga_PAFTOL025241",
                       "rauvolfioids_Tabernaemontaneae_Tabernaemontaninae_Voacanga_thouarsii_PAFTOL006637",
                       "rauvolfioids_Vinceae_Catharanthinae_Kamettia_chandeei_PAFTOL025465",
                       "rauvolfioids_Vinceae_Kopsiinae_Kopsia_arborea_GAP022171",
                       "rauvolfioids_Vinceae_Rauvolfiinae_Rauvolfia_vomitoria_PAFTOL006556",
                       "rauvolfioids_Vinceae_Vincinae_Vinca_difformis_PAFTOL004938",
                       "rauvolfioids_Willughbeieae_Lacmelleinae_Couma_rigida_PAFTOL025227",
                       "rauvolfioids_Willughbeieae_Lacmelleinae_Hancornia_speciosa_PAFTOL025229",
                       "rauvolfioids_Willughbeieae_Lacmelleinae_Lacmellea_bahiensis_PAFTOL025231",
                       "rauvolfioids_Willughbeieae_Lacmelleinae_Parahancornia_negroensis_PAFTOL025499",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Ancylobothrys_amoena_PAFTOL031719",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Chamaeclitandra_henriquesiana_PAFTOL025175",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Clitandra_cymulosa_PAFTOL025475",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Cylindropsis_parvifolia_PAFTOL025307",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Dictyophleba_lucida_PAFTOL025391",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Landolphia_incerta_PAFTOL004111",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Orthopichonia_barteri_PAFTOL031722",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Orthopichonia_seretii_PAFTOL025253",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Pacouria_guianensis_PAFTOL025233",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Saba_comorensis_PAFTOL025281",
                       "rauvolfioids_Willughbeieae_Landolphiinae_Vahadenia_caillei_PAFTOL025357",
                       "rauvolfioids_Willughbeieae_Leuconotidinae_Bousigonia_mekongensis_PAFTOL025503",
                       "rauvolfioids_Willughbeieae_Leuconotidinae_Cyclocotyla_congolensis_PAFTOL025255",
                       "rauvolfioids_Willughbeieae_Leuconotidinae_Leuconotis_anceps_PAFTOL025339",
                       "rauvolfioids_Willughbeieae_Willughbeiinae_Willughbeia_angustifolia_PAFTOL025235")

# Define apocynoids
apocynoids_taxa <- c("apocynoids_Apocyneae_Amphineuriinae_Sindechites_henryi_PAFTOL031695",
                     "apocynoids_Apocyneae_Apocyninae_Apocynum_venetum_PAFTOL004112",
                     "apocynoids_Apocyneae_Beaumontiinae_Vallaris_solanacea_PAFTOL004042",
                     "apocynoids_Apocyneae_Chonemorphinae_Amalocalyx_microlobus_PAFTOL031697",
                     "apocynoids_Apocyneae_Chonemorphinae_Chonemorpha_fragrans_PAFTOL031730",
                     "apocynoids_Apocyneae_Chonemorphinae_Micrechites_rhombifolius_PAFTOL031717",
                     "apocynoids_Apocyneae_Ichnocarpinae_Aganosma_heynei_PAFTOL031726",
                     "apocynoids_Apocyneae_Ichnocarpinae_Baharuia_gracilis_PAFTOL031816",
                     "apocynoids_Apocyneae_Ichnocarpinae_Epigynum_auritum_PAFTOL031703",
                     "apocynoids_Apocyneae_Ichnocarpinae_Pottsia_laxiflora_PAFTOL031715",
                     "apocynoids_Apocyneae_Papuechitinae_Anodendron_affine_PAFTOL031727",
                     "apocynoids_Apocyneae_Papuechitinae_Papuechites_aambe_PAFTOL025337",
                     "apocynoids_Baisseeae_nosubtribe_Baissea_sp_PAFTOL004097",
                     "apocynoids_Baisseeae_nosubtribe_Dewevrella_cochliostema_PAFTOL031758",
                     "apocynoids_Baisseeae_nosubtribe_Motandra_paniculata_PAFTOL025239",
                     "apocynoids_Baisseeae_nosubtribe_Oncinotis_tenuiloba_PAFTOL031812",
                     "apocynoids_EOM_Echitinae_Asketanthera_picardae_PAFTOL025329",
                     "apocynoids_EOM_Echitinae_Echites_umbellatus_PAFTOL004048",
                     "apocynoids_EOM_Echitinae_Thenardia_galeottiana_PAFTOL031814",
                     "apocynoids_EOM_Laubertiinae_Hylaea_arborescens_PAFTOL025493",
                     "apocynoids_EOM_Laubertiinae_Laubertia_contorta_PAFTOL031706",
                     "apocynoids_EOM_nosubtribe_Cycladenia_humilis_PAFTOL031773",
                     "apocynoids_EOM_nosubtribe_Elytropus_chilensis_PAFTOL031763",
                     "apocynoids_EOM_nosubtribe_Forsteronia_acouci_PAFTOL004047",
                     "apocynoids_EOM_nosubtribe_Pinochia_corymbosa_PAFTOL031708",
                     "apocynoids_EOM_nosubtribe_Stipecoma_peltigera_PAFTOL031701",
                     "apocynoids_EOM_Parsonsiinae_Artia_balansae_PAFTOL025405",
                     "apocynoids_EOM_Peltastinae_Temnadenia_odorifera_PAFTOL031711",
                     "apocynoids_EOM_Pentalinoninae_Angadenia_berteroi_PAFTOL031748",
                     "apocynoids_EOM_Pentalinoninae_Pentalinon_luteum_PAFTOL031778",
                     "apocynoids_Malouetieae_Galactophorinae_Galactophora_schomburgkiana_PAFTOL031704",
                     "apocynoids_Malouetieae_Malouetiinae_Allowoodsonia_whitmorei_PAFTOL025437",
                     "apocynoids_Malouetieae_Malouetiinae_Eucorymbia_alba_PAFTOL031716",
                     "apocynoids_Malouetieae_Malouetiinae_Funtumia_africana_PAFTOL031735",
                     "apocynoids_Malouetieae_Malouetiinae_Kibatalia_macrophylla_PAFTOL031733",
                     "apocynoids_Malouetieae_Malouetiinae_Mascarenhasia_lisianthiflora_PAFTOL005311",
                     "apocynoids_Malouetieae_Malouetiinae_Spirolobium_cambodianum_PAFTOL031710",
                     "apocynoids_Nerieae_Alafiinae_Farquharia_elliptica_PAFTOL031691",
                     "apocynoids_Nerieae_Alafiinae_Isonema_smeathmannii_PAFTOL031767",
                     "apocynoids_Nerieae_Alafiinae_Strophanthus_eminii_PAFTOL031743",
                     "apocynoids_Nerieae_Neriinae_Adenium_obesum_PAFTOL004037",
                     "apocynoids_Nerieae_Neriinae_Nerium_oleander_PAFTOL031713",
                     "apocynoids_Rhabdadenieae_nosubtribe_Rhabdadenia_madida_PAFTOL005330",
                     "apocynoids_Wrightieae_nosubtribe_Pleioceras_barteri_PAFTOL031782",
                     "apocynoids_Wrightieae_nosubtribe_Stephanostema_stenocarpum_PAFTOL005310")

# Define Periplocoideae
Periplocoideae_taxa <- c("Periplocoideae_notribe_nosubtribe_Atherandra_acutifolia_PAFTOL031809",
                         "Periplocoideae_notribe_nosubtribe_Baroniella_camptocarpoides_PAFTOL031749",
                         "Periplocoideae_notribe_nosubtribe_Baseonema_gregorii_PAFTOL031750",
                         "Periplocoideae_notribe_nosubtribe_Batesanthus_purpureus_PAFTOL025247",
                         "Periplocoideae_notribe_nosubtribe_Buckollia_volubilis_PAFTOL031751",
                         "Periplocoideae_notribe_nosubtribe_Camptocarpus_mauritianus_PAFTOL031752",
                         "Periplocoideae_notribe_nosubtribe_Chlorocyathus_monteiroae_PAFTOL025505",
                         "Periplocoideae_notribe_nosubtribe_Cryptolepis_nigrescens_PAFTOL025303",
                         "Periplocoideae_notribe_nosubtribe_Cryptostegia_madagascariensis_PAFTOL025359",
                         "Periplocoideae_notribe_nosubtribe_Ectadium_virgatum_PAFTOL031762",
                         "Periplocoideae_notribe_nosubtribe_Epistemma_assianum_PAFTOL031765",
                         "Periplocoideae_notribe_nosubtribe_Finlaysonia_obovata_GAP022109",
                         "Periplocoideae_notribe_nosubtribe_Ischnolepis_natalensis_PAFTOL027331",
                         "Periplocoideae_notribe_nosubtribe_Maclaudia_felixii_PAFTOL031820",
                         "Periplocoideae_notribe_nosubtribe_Mondia_whitei_PAFTOL031737",
                         "Periplocoideae_notribe_nosubtribe_Myriopteron_extensum_PAFTOL025205",
                         "Periplocoideae_notribe_nosubtribe_Pentopetia_albicans_PAFTOL031740",
                         "Periplocoideae_notribe_nosubtribe_Phyllanthera_grayi_GAP022225",
                         "Periplocoideae_notribe_nosubtribe_Phyllanthera_nymanii_PAFTOL025419",
                         "Periplocoideae_notribe_nosubtribe_Raphionacme_hirsuta_PAFTOL031783",
                         "Periplocoideae_notribe_nosubtribe_Sacleuxia_newii_PAFTOL031786",
                         "Periplocoideae_notribe_nosubtribe_Sarcorrhiza_epiphytica_PAFTOL025473",
                         "Periplocoideae_notribe_nosubtribe_Schlechterella_abyssinica_PAFTOL027333",
                         "Periplocoideae_notribe_nosubtribe_Stomatostemma_monteiroae_PAFTOL027339",
                         "Periplocoideae_notribe_nosubtribe_Streptocaulon_juventas_PAFTOL031696",
                         "Periplocoideae_notribe_nosubtribe_Tacazzea_venosa_PAFTOL004103")

# Define Secamonoideae
Secamonoideae_taxa <- c("Secamonoideae_notribe_nosubtribe_Calyptranthera_gautieri_PAFTOL031804",
                        "Secamonoideae_notribe_nosubtribe_Genianthus_micranthus_PAFTOL031821",
                        "Secamonoideae_notribe_nosubtribe_Pervillaea_tomentosa_PAFTOL031707",
                        "Secamonoideae_notribe_nosubtribe_Secamone_schweinfurthii_PAFTOL005226",
                        "Secamonoideae_notribe_nosubtribe_Secamonopsis_madagascariensis_PAFTOL025361",
                        "Secamonoideae_notribe_nosubtribe_Toxocarpus_kleinii_PAFTOL031792",
                        "Secamonoideae_notribe_nosubtribe_Trichosandra_borbonica_PAFTOL031803")

# Define Asclepiadoideae
Asclepiadoideae_taxa <- c("Asclepiadoideae_Asclepiadeae_Asclepiadinae_Asclepias_barjoniifolia_PAFTOL004064",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Aspidoglossum_biflorum_PAFTOL031771",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Aspidonepsis_cognata_PAFTOL031798",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Calciphila_gillettii_PAFTOL027235",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Cordylogyne_globosa_PAFTOL027237",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Fanninia_caloglossa_PAFTOL027239",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Glossostelma_spathulatum_PAFTOL031799",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Kanahia_laniflora_PAFTOL031769",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Margaretta_rosea_PAFTOL031807",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Miraglossum_pulchellum_PAFTOL031898",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Pachycarpus_grandiflorus_PAFTOL031775",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Parapodium_costatum_PAFTOL031777",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Periglossum_mkenii_PAFTOL031693",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Schizoglossum_atropurpureum_PAFTOL031806",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Solenostemma_oleifolium_PAFTOL031742",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Stathmostelma_gigantiflorum_PAFTOL031802",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Stenostelma_capense_PAFTOL031790",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Woodia_mucronata_PAFTOL027241",  
                          "Asclepiadoideae_Asclepiadeae_Asclepiadinae_Xysmalobium_undulatum_PAFTOL031746",  
                          "Asclepiadoideae_Asclepiadeae_Astephaninae_Microloma_sagittatum_PAFTOL031770",  
                          "Asclepiadoideae_Asclepiadeae_Astephaninae_Oncinema_lineare_PAFTOL031774",  
                          "Asclepiadoideae_Asclepiadeae_Cynanchinae_Cynanchum_gracillimum_PAFTOL025271",  
                          "Asclepiadoideae_Asclepiadeae_Cynanchinae_Decanema_bojerianum_PAFTOL031789",  
                          "Asclepiadoideae_Asclepiadeae_Cynanchinae_Schizostephanus_alatus_PAFTOL031787",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Anemotrochus_eggersii_PAFTOL031699",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Chloropetalum_denticulatum_PAFTOL031729",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Cristobalia_bella_PAFTOL027243",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Fischeria_scandens_PAFTOL031817",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Gonolobus_barbatus_PAFTOL031800",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Ibatia_demuneri_PAFTOL027245",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Lachnostoma_ecuadorense_PAFTOL027247",
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Macroscepis_sp_PAFTOL031738",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Peruviasclepias_aliciae_PAFTOL031791",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Phaeostemma_kelleri_PAFTOL031813",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Pherotrichis_sp_PAFTOL031757",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Polystemma_scopulorum_PAFTOL031694",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Pruskortizia_macrocarpa_PAFTOL027249",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Rhytidostemma_sp_PAFTOL027335",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Riparoampelos_amazonicus_PAFTOL031818",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Rojasia_gracilis_PAFTOL027337",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Suberogerens_cyclophylla_PAFTOL031822",  
                          "Asclepiadoideae_Asclepiadeae_Gonolobinae_Tylodontia_stipitata_PAFTOL031714",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Barjonia_laxa_CB448",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Blepharodon_ampliflorum_CB9",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Blepharodon_lineare_CB8",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Ditassa_auriflora_CB95",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Ditassa_blanchetti_CB92",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Ditassa_grandiflora_PAFTOL025209",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Ditassa_mucronata_CB90",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Ditassa_obcordata_CB89",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Ditassa_taxifolia_CB93",  
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Hemipogon_abietoides_CB81",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Hemipogon_acerosus_CB453",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Metastelma_ditassoides_CB82",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Metastelma_harley_CB10",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Metastelma_parviflorum_CB83",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_acerosa_CB411",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_cordata_CB198",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_decussata_CB205",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_ditassoides_CB207",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_grazielae_CB449",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_harley_CB206",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_hemipogonoides_CB208",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_lourtegiae_CB433",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_magisteriana_CB435",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_micromeria_CB436",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_parva_CB443",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_refractifolia_CB442",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_semirii_CB444",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Minaria_volubilis_CB216",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Morilloa_lutea_PAFTOL025211",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Nephradenia_acerosa_CB451",
                          "Asclepiadoideae_Asclepiadeae_Metastelmatinae_Nephradenia_asparagoides_CB452",
                          "Asclepiadoideae_Asclepiadeae_Orthosiinae_Jobinia_connivens_PAFTOL025213",
                          "Asclepiadoideae_Asclepiadeae_Orthosiinae_Monsanima_morrenioides_PAFTOL025215",
                          "Asclepiadoideae_Asclepiadeae_Orthosiinae_Orthosia_dusenii_PAFTOL025217",
                          "Asclepiadoideae_Asclepiadeae_Orthosiinae_Orthosia_scoparia_PAFTOL025207",
                          "Asclepiadoideae_Asclepiadeae_Orthosiinae_Scyphostelma_harlingii_PAFTOL031709",
                          "Asclepiadoideae_Asclepiadeae_Oxypetalinae_Funastrum_angustissimum_PAFTOL031776",
                          "Asclepiadoideae_Asclepiadeae_Oxypetalinae_Philibertia_solanoides_PAFTOL031780",
                          "Asclepiadoideae_Asclepiadeae_Oxypetalinae_Tweedia_brunonis_PAFTOL027251",
                          "Asclepiadoideae_Asclepiadeae_Pentacyphinae_Pentacyphus_camargoi_PAFTOL027253",
                          "Asclepiadoideae_Asclepiadeae_Tassadiinae_Tassadia_stricta_CB457",
                          "Asclepiadoideae_Asclepiadeae_Topeinae_Topea_micrantha_PAFTOL027255",
                          "Asclepiadoideae_Asclepiadeae_Tylophorinae_Pentatropis_capensis_PAFTOL031739",
                          "Asclepiadoideae_Ceropegieae_Anisotominae_Anisotoma_cordifolia_PAFTOL031772",
                          "Asclepiadoideae_Ceropegieae_Anisotominae_Emplectanthus_cordatus_PAFTOL031702",
                          "Asclepiadoideae_Ceropegieae_Anisotominae_Neoschumannia_kamerunensis_PAFTOL031897",
                          "Asclepiadoideae_Ceropegieae_Anisotominae_Riocreuxia_torulosa_PAFTOL031785",
                          "Asclepiadoideae_Ceropegieae_Anisotominae_Sisyranthus_barbatus_PAFTOL025439",
                          "Asclepiadoideae_Ceropegieae_Heterostemminae_Heterostemma_herbertii_PAFTOL031705",
                          "Asclepiadoideae_Ceropegieae_Leptadeniinae_Conomitra_linearis_PAFTOL031808",
                          "Asclepiadoideae_Ceropegieae_Leptadeniinae_Leptadenia_madagascariensis_PAFTOL031797",
                          "Asclepiadoideae_Ceropegieae_Leptadeniinae_Orthanthera_jasminiflora_PAFTOL027257",
                          "Asclepiadoideae_Ceropegieae_Leptadeniinae_Pentasacme_caudatum_PAFTOL027259",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Apteranthes_tuberculata_PAFTOL027263",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Australluma_peschii_PAFTOL027265",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Baynesia_lophophora_PAFTOL027267",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Caralluma_dalzielii_PAFTOL031824",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Ceropegia_sandersonii_PAFTOL005222",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Desmidorchis_retrospiciens_PAFTOL027273",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Duvalia_elegans_PAFTOL031760",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Duvaliandra_dioscoridis_PAFTOL027275",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Echidnopsis_cereiformis_PAFTOL031761",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Edithcolea_grandis_PAFTOL027277",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Huernia_barbata_PAFTOL031766",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Larryleachia_perlata_PAFTOL027279",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Lavrania_haagnerae_PAFTOL027281",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Monolluma_socotrana_PAFTOL027283",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Ophionella_arcuata_PAFTOL027285",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Orbea_hardyi_PAFTOL027287",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Pectinaria_articulata_PAFTOL027289",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Piaranthus_geminatus_PAFTOL027291",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Pseudolithos_dodsonianus_PAFTOL027293",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Pseudolithos_mccoyi_PAFTOL027261",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Quaqua_parviflora_PAFTOL031801",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Rhytidocaulon_macrolobum_PAFTOL027295",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Richtersveldia_columnaris_PAFTOL027297",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Socotrella_dolichocnema_PAFTOL027299",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Stapelianthus_decaryi_PAFTOL027301",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Stapeliopsis_saxatilis_PAFTOL027303",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Tavaresia_barklyi_PAFTOL027305",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Tridentea_gemmiflora_PAFTOL027307",
                          "Asclepiadoideae_Ceropegieae_Stapeliinae_Tromotriche_ruschiana_PAFTOL027309",
                          "Asclepiadoideae_Eustegieae_nosubtribe_Emicocarpus_fissifolius_PAFTOL025387",
                          "Asclepiadoideae_Eustegieae_nosubtribe_Eustegia_minuta_PAFTOL004041",
                          "Asclepiadoideae_Fockeeae_nosubtribe_Cibirhiza_dhofarensis_PAFTOL031755",
                          "Asclepiadoideae_Fockeeae_nosubtribe_Fockea_edulis_PAFTOL005223",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Anisopus_mannii_PAFTOL031795",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Campestigma_purpureum_PAFTOL031700",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Cionura_erecta_PAFTOL031756",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Cosmostigma_cordatum_PAFTOL025333",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Dalzielia_oblanceolata_PAFTOL025355",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Dischidanthus_urceolatus_PAFTOL027311",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Dischidia_nummularia_GAP022127",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Dregea_floribunda_PAFTOL031759",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Gongronema_filipes_PAFTOL027313",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Gymnemopsis_calcicola_PAFTOL025463",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Harmandiella_cordifolia_PAFTOL027317",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Jasminanthes_suaveolens_PAFTOL031768",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Lygisma_inflexum_PAFTOL027321",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Marsdenia_angolensis_PAFTOL027315",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Oreosparte_celebica_PAFTOL027323",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Papuahoya_bykulleana_PAFTOL027325",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Pycnorhachis_maingayi_PAFTOL025377",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Rhyssolobium_dumosum_PAFTOL031784",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Ruehssia_carvalhoi_PAFTOL004072",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Sarcolobus_vittatus_GAP022293",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Sicyocarpus_verrucosus_PAFTOL027329",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Stigmatorhynchus_umbellifer_PAFTOL031698",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Treutlera_insignis_PAFTOL025399",
                          "Asclepiadoideae_Marsdenieae_nosubtribe_Wattakaka_volubilis_PAFTOL031745")


# --- Utility functions ---
assign_tip_colors <- function(tips) {
  pal <- as.character(paletteer::paletteer_d("rcartocolor::Bold")[1:7])  # Use only the 7 defined
  vec <- rep(NA, length(tips))
  names(vec) <- tips
  vec[tips %in% rogue_taxa] <- "red"
  vec[tips %in% outgroups_taxa] <- "darkgray"
  vec[tips %in% rauvolfioids_taxa] <- pal[1]
  vec[tips %in% apocynoids_taxa] <- pal[2]
  vec[tips %in% Periplocoideae_taxa] <- pal[3]
  vec[tips %in% Secamonoideae_taxa] <- pal[4]
  vec[tips %in% Asclepiadoideae_taxa] <- pal[5]
  return(vec)
}


PP_calc <- function(AS_A) {
  PP1 <- as.numeric(AS_A@data$pp1) * 100
  PP <- as.data.frame(PP1)
  PP$PP2 <- as.numeric(AS_A@data$pp2) * 100
  PP$PP3 <- as.numeric(AS_A@data$pp3) * 100
  PP$node <- AS_A@data$node
  return(PP)
}

plot_tree_with_ppl <- function(AS_tree, title) {
  phylo <- AS_tree@phylo
  tip_colors <- assign_tip_colors(phylo$tip.label)
  
  base <- ggtree(phylo, ladderize = TRUE, branch.length = "none") +
    geom_tree(color = "darkgray") +  # All branches black
    geom_tiplab(aes(label = label, color = label), size = 3.0, hjust = 0) +
    xlim_tree(70) +
    scale_color_manual(values = tip_colors, guide = "none") +
    ggtitle(title)
  
  pp_support <- PP_calc(AS_tree)
  pies <- nodepie(pp_support, cols = 1:3, color = c(PP1 = '#006CD1', PP2 = '#994F00', PP3 = 'gray'))
  return(inset(base, pies, width = 0.06, height = 0.06))
}

# --- Load data and trees ---
data <- read.csv("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Data/names_dt1.csv")

rename_tree <- function(tree_path) {
  tr <- read.astral(tree_path)
  rename_taxa(tr, data, key = 1, value = 2)
}

AS_all_raxml_ppl <- rename_tree("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Data/dataset2_diamond/all_raxml_AAtoDNA/SpeciesTree_PxrrRooted_raxml_allebg.tre")
AS_nha_raxml_ppl <- rename_tree("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Data/dataset2_diamond/nha_raxml_AAtoDNA/SpeciesTree_PxrrRooted_raxml_nohard.tre")

# --- Plot trees with Q pie charts and tip coloring ---
p_all_raxml_ppl <- plot_tree_with_ppl(AS_all_raxml_ppl, "All genes RAxML (n=332)")
p_nha_raxml_ppl <- plot_tree_with_ppl(AS_nha_raxml_ppl, "Only good genes RAxML (n=212)")

# --- Export PDF ---
pdf("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Results/Supplement_2_PPLpies.pdf", 20, 50)
grid.arrange(p_all_raxml_ppl, p_nha_raxml_ppl, nrow = 1)
dev.off()
