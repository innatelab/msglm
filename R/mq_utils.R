# MaxQuant data import utils
# 
# Author: astukalov
###############################################################################

source( file.path( base_scripts_path, 'misc/reshape.R' ) )

require(readr)
require(tidyr)
require(stringr)
require(CeMMmisc)
require(quantreg)

recombine_dlms <- function(dlms, sep=";") {
    dlms_expanded <- unique(unlist(lapply(dlms, strsplit, sep, fixed=TRUE)))
    paste0(dlms_expanded[!is.na(dlms_expanded)], collapse=sep)
}

expand_protgroups <- function(protgroup_ids) {
    protgroups2protgroup.list <- strsplit(protgroup_ids, ';', fixed=TRUE)
    data.frame( protgroup_ids = rep.int(protgroup_ids, sapply(protgroups2protgroup.list, length)),
                protgroup_id = as.integer(unlist(protgroups2protgroup.list)),
                stringsAsFactors = FALSE )
}

selectUniprotACs <- function( acs, valid_acs )
{
    acs_noiso <- stripUniprotIsoform( acs )
    if ( !is.null( valid_acs ) ) {
        acs_noiso <- intersect( acs_noiso, valid_acs )
    }
    if ( length(acs_noiso) == 1 ) {
        return ( acs_noiso ) }
    } else if ( length(acs_noiso) > 1 ) {
        # return the first "classical" UniProt AC starting with O, P or Q 
        is_classical_uprot <- str_detect(acs_noiso, '^[POQ]')
        return ( acs_noiso[[ order(!is_classical_uprot, acs_noiso)[[1]] ]] )
    } else { return ( NA_character_ ) }
}

read.MaxQuant <- function( filename, layout = c( "wide", "long" ),
        row_id_cols = c( 'protein_ac_noiso' ),
        protein_ac_cols, protein_info = NULL,
        measures.regex )
{
      res <- read.table( filename, sep = '\t', header = TRUE, check.names = FALSE, stringsAsFactors = FALSE )
      #print(colnames(res))
      measures <- gsub( '\\\\', '', measures.regex )
      measures.fixed <- measures %>% gsub( '\\s*\\[\\%\\]', '', . ) %>% gsub( '/', '', . ) %>% gsub( '\\+', 'and', . ) %>% gsub( '\\s', '_', . )
      names(measures.fixed) <- measures

      # separate measure from SILAC label and from the run label and fix the total/summary names 
      colnames( res ) <- colnames( res ) %>%
                sub( '^Sequence coverage (.+) \\[\\%\\]$', 'Sequence coverage [%] \\1', . ) %>% # fix the position of [%] to make the naming scheme consitent 
                sub( paste0('^(',paste0( measures.regex, collapse='|'),')(\\s([HLM]))?(\\s([^HLM].*))?$' ), '\\1.\\3__\\5', . ) %>%
                sub( '\\.__', '.Sum__', . ) %>% sub( '\\__$', '\\__Everything', . )
      #print(colnames(res))

      # fix UniProt ACs
      for ( i in seq_along( protein_ac_cols ) ) {
            if ( !( protein_ac_cols[[i]] %in% colnames(res) ) ) {
                  warning( 'Protein AC column `', protein_ac_cols[[i]], '` not found' )
                  next
            }
            fixed_acs <- sapply( strsplit( res[[ protein_ac_cols[[i]] ]], ';', fixed = TRUE ), selectUniprotACs,
                if ( is.null(protein_info) ) NULL else protein_info$protein_ac_noiso )
            na_mask <- is.na(fixed_acs)
            if ( any( na_mask ) ) {
                  warning( sum( na_mask ), ' records have no valid protein AC, removed: ',
                                   paste0( res[[ protein_ac_cols[[i]] ]][na_mask], collapse = ' ' ) )
                  res <- res[ !na_mask, , drop = FALSE ]
            }
            res[,names(protein_ac_cols)[[i]] ] <- fixed_acs[!na_mask] 
      }

      quant.columns.list <- lapply( measures.regex, function( mes ) grep( paste0( '^', mes, '\\.' ), colnames( res ), value = TRUE ) )
      names(quant.columns.list) <- measures
      quant.columns <- unlist( quant.columns.list )
      channels <- sort( unique( gsub( '^[^.]+\\.', '', quant.columns ) ) )
      channels.df <- cbind( mschannel = channels, do.call( rbind, strsplit( channels, '__', fixed = TRUE ) ) ) %>%
                as.data.frame( stringsAsFactors = FALSE )
      colnames(channels.df) <- c( 'mschannel', 'mstag', 'msrun' )
      channels.df <- dplyr::mutate(channels.df,
        mstag = factor(mstag, levels = intersect(c('H', 'M', 'L', 'Sum'), unique(mstag))),
        mstag_type = if_else(mstag %in% c('H', 'M', 'L'), 'SILAC', 'label_free'))

      quant.columns.df <- expand.grid( measure = measures, mschannel = channels,
                    stringsAsFactors = FALSE ) %>%
                mutate( colname = paste0( measure, '.', channel ),
                                measure_fixed = measures.fixed[ measure ],
                                colname_fixed = paste0( measure_fixed, '.', mschannel ) ) %>%
                inner_join( channels.df )

      #print(str(res))
      #print( quant.columns.df )
      res <- res %>% mutate_each_( funs(as.numeric), backquote(quant.columns) ) %>%
                group_by_( .dots = paste0( '`', row_id_cols, '`' ) ) %>%
                summarise_each_( funs(max), backquote(quant.columns) ) %>%
                dplyr::ungroup()

      if ( match.arg( layout ) == 'long' ) {
            quant.columns.df$exists <- quant.columns.df$colname %in% colnames(res)
            res[, subset(quant.columns.df,!exists)$colname ] <- NA # add missing column names to make it reshapeable
            rename_cols <- paste0( '`', quant.columns.df$colname, '`' )
            names(rename_cols) <- quant.columns.df$colname_fixed
            res <- rename_( res, .dots = rename_cols )
            quant.columns.full_list <- lapply( measures.fixed, function(mes) subset(quant.columns.df,measure_fixed==mes)$colname_fixed )
            names(quant.columns.full_list) <- measures.fixed
            res <- reshape( res, direction = 'long',
                idvar = row_id_cols,
                varying = quant.columns.full_list,
                v.names = names(quant.columns.full_list),
                timevar = 'mschannel',
                times = sort( unique( quant.columns.df$mschannel ) ),
                sep = '.' )
            res <- inner_join( res, channels.df ) %>% mutate( mschannel = NULL ) # add mstags and msrun info
            # fix NA
            for ( mes in measures.fixed ) {
                  res[ !is.na( res[[mes]] ) & res[[mes]] == 0, mes ] <- NA
            }
      }
      if ( !is.null( protein_info ) ) {
            prot_info_cols <- c('protein_label','genename','molweight','contaminant_family')
            join_by <- c( protein_ac_noiso = 'protein_ac_noiso' )
            names(join_by) <- names(protein_ac_cols)[[1]]
            res <- res %>% left_join( dplyr::select_( protein_info, .dots = c( 'protein_ac_noiso', prot_info_cols ) ),
                by = join_by )
      }
      return ( res )
}

gsub_columns <- function(df, from, to) {
  colnames(df) <- gsub(from, to, colnames(df))
  df
}

read.MaxQuant.ProteinGroups <- function( folder_path, file_name = 'proteinGroups.txt',
                                         protein_info = NULL, layout = c( "wide", "long" ),
                                         nrows = Inf, import_data = c(), guess_max = min(10000L, nrows) )
{
    proteinGroups.df <- read_tsv( file.path( folder_path, file_name ), n_max = nrows,
                             col_types = cols( `Fasta headers` = "c" ),
                             na = c( "", "NA", "n. def.", "n.def." ), guess_max = guess_max )
    col_renames <- c("protgroup_id" = "id", "protein_acs" = "Protein IDs",
                     "majority_protein_acs" = "Majority protein IDs",
                     "gene_names" = "Gene names", "protein_names" = "Protein names",
                     "fasta_headers" = "Fasta headers", "n_proteins" = "Number of proteins",
                     "score" = "Score", "q_value" = "Q-value",
                     "seqlen" = "Sequence length", "seqlens" = "Sequence lengths",
                     "mol_weight_kDa" = "Mol. weight [kDa]", "seqcov" = "Sequence coverage [%]",
                     "unique_razor_seqcov" = "Unique + razor sequence coverage [%]",
                     "is_contaminant" = "Potential contaminant", "is_reverse" = "Reverse")
    res.df <- dplyr::select_(proteinGroups.df, .dots = backquote(col_renames)[col_renames %in% colnames(proteinGroups.df)]) %>%
      dplyr::mutate( is_contaminant = !is.na(is_contaminant) & is_contaminant == '+',
                     is_reverse = !is.na(is_reverse) & is_reverse == '+' )
    col_info <- list( protgroup = colnames(res.df) )
    if ('intensity' %in% import_data) {
      intensities.df <- proteinGroups.df %>% dplyr::select(starts_with("Intensity")) %>%
        gsub_columns("^Intensity\\s([LMH](\\s|$))", "Intensity.\\1") %>%
        gsub_columns("^Intensity(\\s|$)", "Intensity.Sum\\1") %>%
        mutate_each(., funs(ifelse(.==0.0, NA, .)))
      res.df <- bind_cols(res.df, intensities.df)
      col_info$intensity <- colnames(intensities.df)
    }
    if ('LFQ' %in% import_data) {
      lfq.df <- proteinGroups.df %>% dplyr::select(starts_with("LFQ intensity")) %>%
         gsub_columns("^LFQ intensity\\s([LMH](\\s|$))", "LFQ_Intensity.\\1") %>%
         gsub_columns("^LFQ intensity(\\s|$)", "LFQ_Intensity.Sum\\1") %>%
         mutate_each(., funs(ifelse(.==0.0, NA, .)))
      res.df <- bind_cols(res.df, lfq.df)
      col_info$LFQ <- colnames(lfq.df)
    }
    if ('iBAQ' %in% import_data) {
      ibaq.df <- proteinGroups.df %>% dplyr::select(starts_with("iBAQ")) %>%
        gsub_columns("^iBAQ\\s([LMH](\\s|$))", "iBAQ.\\1") %>%
        gsub_columns("^iBAQ(\\s|$)", "iBAQ.Sum\\1") %>%
        mutate_each(., funs(ifelse(.==0.0, NA, .)))
      res.df <- bind_cols(res.df, ibaq.df)
      col_info$iBAQ <- colnames(ibaq.df)
    }
    if ('ratio' %in% import_data) {
      ratios.df <- proteinGroups.df %>% dplyr::select(starts_with("Ratio")) %>%
        gsub_columns("^Ratio\\s", "Ratio.")
      res.df <- bind_cols(res.df, ratios.df)
      col_info$ratio <- colnames(ratios.df)
    }
    attr(res.df, "column_groups") <- col_info
    return (res.df)
}

read.MaxQuant.Peptides <- function( folder_path, file_name = 'peptides.txt',
                                    import_data = c(), nrows = Inf )
{
  peptides.df <- read_tsv( file.path( folder_path, file_name ), n_max = nrows,
                           col_types = cols( `Proteins` = "c",
                                             `Protein group IDs` = "c",
                                             `Mod. peptide IDs` = "c",
                                             `Leading razor protein` = "c",
                                             `Charges` = "c" ),
                           na = c( "", "NA", "n. def.", "n.def." ),
                           guess_max = 20000L ) %>%
    dplyr::mutate( peptide_id = row_number()-1,
                   is_reverse = !is.na(`Reverse`) & `Reverse` == '+',
                   is_contaminant = !is.na(`Potential contaminant`) & `Potential contaminant` == '+',
                   is_shared_by_groups = is.na(`Unique (Groups)`) | !(`Unique (Groups)` %in% c('yes', '+')),
                   is_shared = is_shared_by_groups,
                   is_shared_by_proteins = is.na(`Unique (Proteins)`) | !(`Unique (Proteins)` %in% c('yes', '+')) )
  #message(paste0(colnames(peptides.df), collapse='\n'))
  res.df <- peptides.df %>%
    dplyr::select(peptide_id, seq = Sequence, seq_len = Length, n_miscleavages = `Missed cleavages`,
                  nterm_window = `N-term cleavage window`, cterm_window = `C-term cleavage window`,
                  aa_before = `Amino acid before`, aa_after = `Amino acid after`,
                  aa_first = `First amino acid`, aa_last = `Last amino acid`,
                  protgroup_ids = `Protein group IDs`,
                  pepmod_ids = `Mod. peptide IDs`,
                  protein_acs = `Proteins`, lead_razor_protein_ac = `Leading razor protein`,
                  start_pos = `Start position`, end_pos = `End position`,
                  mass = Mass, charges = `Charges`, score = Score, PEP = PEP,
                  is_reverse, is_contaminant, is_shared_by_groups, is_shared, is_shared_by_proteins)
  col_info <- list( peptide = colnames(res.df) )
  if ('intensity' %in% import_data) {
    intensities.df <- peptides.df %>% dplyr::select(starts_with("Intensity")) %>%
      gsub_columns("^Intensity\\s([LMH](\\s|$))", "Intensity.\\1") %>%
      gsub_columns("^Intensity(\\s|$)", "Intensity.Sum\\1") %>%
      mutate_each(., funs(ifelse(.==0.0, NA, .)))
    res.df <- bind_cols(res.df, intensities.df)
    col_info$intensity <- colnames(intensities.df)
  }
  return (res.df)
}

# read allPeptides.txt table
read.MaxQuant.AllPeptides <- function( folder_path, file_name = 'allPeptides.txt', nrows = Inf )
{
    allPeptides.df <- read_tsv( file.path( folder_path, file_name ), n_max = nrows,
                           col_types = cols( `Type` = "c",
                                             `Raw file` = "c",
                                             `Resolution` = "n"
                                             ),
                           na = c( "", "NA", "n. def.", "n.def." ),
                           guess_max = 20000L )
    allPeptides.df <- dplyr::rename(allPeptides.df,
                                    pep_type = Type,
                                    raw_file = `Raw file`,
                                    #protgroup_ids = `Protein group IDs`,
                                    #fasta_headers = `Fasta headers`,
                                    dp_modif = `DP Modification`,
                                    dp_mass_diff = `DP Mass Difference`) %>%
      dplyr::mutate(pep_type = factor(pep_type),
                    raw_file = factor(raw_file),
                    #protgroup_ids = factor(protgroup_ids),
                    #fasta_headers = factor(fasta_headers),
                    is_dp = !is.na(dp_mass_diff),
                    dp_modif = factor(dp_modif),
                    dp_mass_delta = as.integer(dp_mass_diff*100)/100)
    # remove less useful memory-hungry columns
    #allPeptides.df$Intensities <- NULL
    #allPeptides.df$dp_mass_delta <- as.integer(allPeptides.df$`DP Mass Difference`*100)/100
    return(allPeptides.df)
}

read.MaxQuant.Sites <- function( folder_path, file_name, nrows = Inf, modif = "Phospho (STY)",
                                 import_data = c() )
{
  data.df <- read_tsv( file.path( folder_path, file_name ),
                        col_names = TRUE, n_max = nrows,
                        col_types = cols(`Protein group IDs` = col_character(), .default = col_guess() ),
                        na = c( "", "NA", "n. def.", "n.def." ), guess_max = 20000L )
  sites.df <- data.df %>%
    dplyr::mutate( is_contaminant = !is.na(`Potential contaminant`) & `Potential contaminant` == '+',
                   is_reverse = !is.na(`Reverse`) & `Reverse` == '+' ) %>%
    dplyr::select(site_id = id,
                  protgroup_ids = `Protein group IDs`,
                  leading_protein_acs = `Leading proteins`,
                  protein_ac = `Protein`,
                  protein_acs = `Proteins`,
                  loc_prob = `Localization prob`,
                  delta_score = `Delta score`,
                  score = `Score`,
                  loc_score = `Score for localization`,
                  peptide_ids = `Peptide IDs`,
                  pepmod_ids = `Mod. peptide IDs`,
                  positions = `Positions`,
                  position = `Position`,
                  positions_in_proteins = `Positions within proteins`,
                  aa = `Amino acid`,
                  seq_context = `Sequence window`,
                  modif_context = `Modification window`,
                  gene_names = `Gene names`,
                  protein_names = `Protein names`,
                  fasta_headers = `Fasta headers`,
                  charge = `Charge`,
                  mass_error_ppm = `Mass error [ppm]`,
                  is_contaminant, is_reverse)
  sites.df["n_of_modifs"] <- data.df[[paste0("Number of ", modif)]]
  col_info <- list( site = colnames(sites.df) )
  res.df <- sites.df
  if ('intensity' %in% import_data) {
    intensities.df <- data.df %>% dplyr::select(starts_with("Intensity")) %>%
      gsub_columns("^Intensity\\s([LMH](\\s|_|$))", "Intensity.\\1") %>%
      gsub_columns("^Intensity(\\s|_|$)", "Intensity.Sum\\1") %>%
      gsub_columns("^Intensity.Sum(.+)(___\\d+)$", "Intensity.Sum\\2\\1") %>%
      mutate_each(., funs(as.numeric)) %>%
      mutate_each(., funs(ifelse(.==0.0, NA_real_, .)))
    res.df <- bind_cols(res.df, intensities.df)
    col_info$intensity <- colnames(intensities.df)
  }
  if ('occupancy' %in% import_data) {
    occupancies.df <- data.df %>% dplyr::select(starts_with("Occupancy")) %>%
      gsub_columns("^(Occupancy\\s|Occupancy ratio|Occupancy error scale\\s)([LMH](\\s|$))", "\\1.\\2") %>%
      gsub_columns("^Occupancy ratio", "occupancy_ratio") %>%
      gsub_columns("^Occupancy error scale", "occupancy_error_scale") %>%
      mutate_each(., funs(as.numeric))
    res.df <- bind_cols(res.df, occupancies.df)
    col_info$occupancy <- colnames(occupancies.df)
  }
  if ('ratio' %in% import_data) {
    ratios.df <- data.df %>% dplyr::select(starts_with("Ratio mode/base")) %>%
      gsub_columns("^Ratio mod/base\\s", "ratio_mod2base.") %>%
      mutate_each(funs(as.numeric))
    res.df <- bind_cols(res.df, ratios.df)
    col_info$ratio <- colnames(ratios.df)
  }
  attr(res.df, "column_groups") <- col_info
  return (res.df)
}

backquote <- function( vars ) {
  res <- paste0( "`", vars, "`" )
  names(res) <- names(vars)
  return (res)
}

read.MaxQuant.Evidence_internal <- function( folder_path, file_name = 'evidence.txt',
                                             nrows = Inf, guess_max = min(20000L, nrows) ) {
  message( 'Reading evidence table...' )
  evidence.df <- read_tsv(file.path( folder_path, file_name ), col_names = TRUE, n_max = nrows,
                          col_types = cols( Resolution = 'n',
                                           `Protein group IDs` = 'c',
                                           `Oxidation (M) site IDs` = 'c',
                                           `Peptide ID` = 'i', `Mod. peptide ID` = 'i', `Charge` = 'i',
                                            Type = 'c', `Raw file` = 'c', Experiment = 'c',
                                            Modifications = 'c', `Labeling State` = 'c',
                                           .default = col_guess() ),
                          na = c( "", "NA", "n. def.", "n.def." ), guess_max = guess_max)

  message( 'Renaming and converting evidence table columns...' )
  # fix column names from different MQ versions
  col_renames <- c( evidence_id = 'id',
                    pepmod_id = 'Mod. peptide ID',
                    raw_file = "Raw file", msrun = "Experiment", ident_type = "Type", label_state = 'Labeling State',
                    mass_error_ppm = "Mass Error [ppm]", mass_error_da = "Mass Error [Da]",
                    mass_error_ppm = "Mass error [ppm]", mass_error_da = "Mass error [Da]",
                    uncalib_mass_error_ppm = "Uncalibrated Mass Error [ppm]",
                    uncalib_mass_error_da = "Uncalibrated Mass Error [Da]",
                    seq = 'Sequence', seq_len = 'Length', count_K = 'K Count', count_R = 'R Count', modifs = 'Modifications',
                    n_miscleaved = 'Missed cleavages',
                    mod_seq = 'Modified sequence', charge = 'Charge',
                    protein_acs = 'Proteins', lead_protein_acs = 'Leading proteins', lead_razor_protein_ac = 'Leading razor protein',
                    gene_names = 'Gene Names', protein_names = 'Protein Names', fasta_headers = 'Fasta headers',
                    is_reverse = 'Reverse',
                    is_contaminant = "Potential contaminant", is_contaminant = 'Contaminant',
                    rt_len = 'Retention length', rt_avg = 'Calibrated retention time', rt_start = 'Calibrated retention time start', rt_finish = 'Calibrated retention time finish',
                    protgroup_ids = 'Protein group IDs', # IDs between different folders do not match
                    # Carbamidomethyl (C) site IDs', Oxidation (M) site IDs
                    peptide_id = 'Peptide ID',
                    mz = 'm/z', ms2_mz = 'MS/MS m/z', mass_da = 'Mass', resolution = 'Resolution',
                    delta_mz_da = "Uncalibrated - Calibrated m/z [Da]",
                    delta_mz_ppm = "Uncalibrated - Calibrated m/z [ppm]"
  )
  if (any(str_detect(colnames(evidence.df), "Intensity "))) {
    # intensity is an aggregation column of labeled intensities
    col_renames = c(col_renames, "Intensity Sum" = "Intensity")
  } else {
    # label-free data, rename column so it gets the tag
    col_renames = c(col_renames, "Intensity F" = "Intensity")
  }
  col_renames <- col_renames[ col_renames %in% colnames(evidence.df) ]
  if ( length( col_renames ) > 0 ) {
    evidence.df <- dplyr::rename_( evidence.df, .dots = backquote(col_renames) )
  }
  evidence.df <- dplyr::mutate(evidence.df,
                               msrun = factor(msrun),
                               raw_file = factor(raw_file),
                               ident_type = factor(ident_type),
                               is_contaminant = !is.na(is_contaminant) & is_contaminant == '+',
                               is_reverse = !is.na(is_reverse) & is_reverse == '+')
  return ( evidence.df )
}

glm_corrected_intensities <- function(intensities.df,
                                      use_mstags=FALSE,
                                      glm_max_factor = 50.0,
                                      glm_reldelta_range = c(-0.9, 0.9)) {
  ndatapoints <- nrow(intensities.df)
  nmsprotocols <- if ('msprotocol' %in% colnames(intensities.df)) n_distinct(intensities.df$msprotocol) else 1L
  msrunXcharge <- dplyr::distinct(dplyr::select(dplyr::filter(intensities.df, observed), msrun, charge))
  # GLM could be run also with ncharges = n_distinct(charge) > 1, but often generates artifacts
  ncharges <- max(table(msrunXcharge$msrun))
  nmsruns <- n_distinct(intensities.df$msrun)
  npepmods <- n_distinct(intensities.df$pepmod_id)
  msrunXtag <- dplyr::distinct(dplyr::select(dplyr::filter(intensities.df, observed), msrun, mstag))
  nmstags <- max(table(msrunXtag$msrun))
  rhs_str <- "0"
  if (nmsruns >= 1) {
    if (nmstags > 1) {
      if (use_mstags) {
        rhs_str <- paste0(rhs_str, " + msrun + msrun:mstag")
      } else { # treat each mschannel independently not considering they are from the same run
        rhs_str <- paste0(rhs_str, " + mschannel")
      }
    } else {
      rhs_str <- paste0(rhs_str, " + msrun")
    }
  }
  if (ncharges > 1) {
    rhs_str <- paste0(rhs_str, if_else(npepmods>1, " + pepmod_id:charge", " + charge"))
  }
  if (npepmods > 1) {
    rhs_str <- paste0(rhs_str, " + pepmod_id")
    if (nmsruns > 1) {
      rhs_str <- paste0(rhs_str, " + msrun:pepmod_id")
    }
  }
  max_intensity_fixed = max(intensities.df$intensity_fixed, na.rm = TRUE)
  intensities.df$intensity_glm <- NA_real_
  intensities.df$glm_rhs <- NA_character_
  if (ndatapoints > 1L && rhs_str != "0 + msrun" && rhs_str != "0 + mschannel") {

      # convert GLM terms into factors to generate proper model and
      # use only the factors present in the given intensities chunk
      intensities_copy.df <- dplyr::mutate(intensities.df,
          intensity_fixed_norm = intensity_fixed/max_intensity_fixed,
          pepmod_id = factor(pepmod_id, ordered = FALSE),
          charge = factor(charge, ordered = FALSE),
          mstag = factor(mstag, ordered = FALSE),
          msprotocol = factor(msprotocol, ordered = FALSE)
      )
      fla <- as.formula(paste0("log(intensity_fixed_norm) ~ ", rhs_str))
      for (ltl in -6:1) { # try different scales
        tryCatch({
          # FIXME GLM assumes one pepmod_id
          glm_res <- rq(fla, data=intensities_copy.df, weights = intensities_copy.df$weight)
          intensities.df$intensity_glm <- exp(predict(glm_res))*max_intensity_fixed
          intensities.df$glm_rhs <- rhs_str
          break
        }, error = function(e) {warning(e)}, finally = invisible() )
      }
  }
  if (!is.na(intensities.df$glm_rhs[1])) {
    # GLM succeeded, now check if its results are valid and correct
    intensities.df <- dplyr::mutate(intensities.df,
                                    intensity_glm_reldelta = (intensity_glm - intensity_fixed) / pmax(intensity_glm, intensity_fixed),
                                    is_valid = (intensity_glm_reldelta >= glm_reldelta_range[[1]]) & (!observed | (intensity_glm_reldelta <= glm_reldelta_range[[2]])))
    if (all(intensities.df$is_valid)) {
      # use GLM predictions correcting for
      # anomalously big GLM predictions (e.g. when charge of higest intensity was missing and got predicted)
      intensities.df <- dplyr::mutate(intensities.df,
                                      intensity_corr = pmin(intensity_glm, glm_max_factor*max_intensity_fixed) )
    } else {
      # GLM generates outliers that do not make sense, use the original intensities
      intensities.df$intensity_corr <-intensities.df$intensity
      intensities.df$rhs_glm <- paste0(intensities.df$rhs_glm, " (invalid)")
    }
    intensities.df <- dplyr::mutate(intensities.df,
                                    intensity_corr_reldelta = (intensity_corr - intensity_fixed) / pmax(intensity_corr, intensity_fixed) )
  } else {
    # GLM failed, use the original intensities
    intensities.df <- dplyr::mutate(intensities.df,
                                    intensity_glm_reldelta = NA_real_,
                                    is_valid = NA,
                                    intensity_corr = intensity,
                                    intensity_corr_reldelta = if_else(!is.na(intensity), 0.0, NA_real_))
  }
  return (intensities.df)
}

process.MaxQuant.Evidence <- function( evidence.df, layout = c( "pepmod_msrun", "pepmod_mschannel" ),
                                       min_pepmod_state_freq = 0.9, min_essential_freq = 0.0,
                                       import_data = c("intensity"),
                                       correct_ratios = TRUE,
                                       correct_by_ratio.ref_label = NA,
                                       mode = c("labeled", "label-free"),
                                       mschannel_annotate.f = NULL,
                                       na_weight = 1E-5, min_intensity = 1E+3,
                                       glm.min_intensity = 50*min_intensity,
                                       glm.max_factor = 50,
                                       glm.reldelta_range = c(-0.99, 0.99) )
{
    quant_mode <- match.arg(mode)
    complete_names <- function( vars ) {
      res <- vars
      res_names <- names(vars)
      res_names[res_names==""] <- vars[res_names==""]
      names(res) <- res_names
      return (res)
    }
    message("Extracting MS runs and MS channels info...")
    ilabels <- factor( c("H", "M", "L", "F", 'Sum'), ordered = TRUE,
                       levels = c("H", "M", "L", "F", 'Sum') )
    intensity_columns.df <- expand.grid( measure = 'intensity', mstag = ilabels ) %>%
      mutate( old_name = paste0('Intensity ', mstag),
              new_name = paste0(measure, '.', mstag),
              type = if_else(mstag == 'Sum', 'aggregate', if_else(old_name %in% colnames(evidence.df), 'measured', 'missing')) ) %>%
      dplyr::filter( type != 'missing' ) %>%
      dplyr::mutate( type = factor(type)) %>%
      dplyr::arrange( mstag )
    ilabels <- intensity_columns.df$mstag # restrict to the labels actually used
    msruns.df <- evidence.df %>% dplyr::select(msrun, raw_file) %>% dplyr::distinct()
    mschannels.df <- expand.grid( raw_file = msruns.df$raw_file,
                                  mstag = ilabels ) %>%
      dplyr::inner_join(msruns.df) %>%
      dplyr::mutate(
          mschannel = interaction(msrun, mstag, drop=TRUE, lex.order=TRUE, sep='_'),
          quant_type = if_else(mstag %in% c('H','M','L'), 'SILAC',
                       if_else(mstag == 'Sum', 'aggregate', 'label_free')))
    if (any(mschannels.df$quant_type == 'SILAC')) {
      message('Evidence table contains SILAC-labeled data')
      mode <- 'labeled'
    } else {
      message('Evidence table contains label-free data')
      model <- 'label-free'
    }
    if (!is.null(mschannel_annotate.f)) {
      # get user-defined mschannel annotations
      annot.df <- mschannel_annotate.f(mschannels.df)
      if (nrow(annot.df) != nrow(mschannels.df)) {
        stop("mschannel annotations has incorrect number of rows")
      }
      mschannels_annot.df <- dplyr::left_join(mschannels.df, annot.df)
      if (nrow(mschannels_annot.df) != nrow(mschannels.df)) {
        stop("mschannel annotations do not correctly match mschannels")
      }
      mschannels.df <- mschannels_annot.df
      if ('is_msrun_used' %in% colnames(mschannels.df)) {
        message('Removing unused user-specified msruns')
        mschannels.df <- dplyr::filter(mschannels.df, !is.na(is_msrun_used) & is_msrun_used) %>%
            # drop unused levels
            dplyr::mutate(msrun = factor(msrun),
                          raw_file = factor(raw_file),
                          mschannel = factor(mschannel))
        evidence.df <- dplyr::filter(evidence.df, as.character(raw_file) %in% mschannels.df$raw_file) %>%
            dplyr::mutate(raw_file = factor(raw_file, levels=levels(mschannels.df$raw_file)),
                          msrun = factor(msrun, levels=levels(mschannels.df$msrun)))
      }
      if (n_distinct(mschannels.df$mschannel) != nrow(mschannels.df)) {
        stop('mschannel ids are not unique')
      }
    }
    mschannels.df <- dplyr::mutate(mschannels.df, quant_type = factor(quant_type))
    # NOTE?: the same mod. peptide Id can have multiple mod. sequences (mod at different poses)
    message('Enumerating pepmod states and summarizing channel intensities...')
    intensities.mtx <- as.matrix(evidence.df[,dplyr::filter(intensity_columns.df, mstag != 'Sum')$old_name])
    intens_cols <- backquote(intensity_columns.df$old_name)
    names(intens_cols) <- intensity_columns.df$new_name
    evidence.df <- mutate(evidence.df, 
                          pepmod_state = interaction(pepmod_id, charge, drop = TRUE, lex.order = TRUE, sep = '_'),
                          `Intensity Sum` = rowSums(intensities.mtx, na.rm = TRUE),
                          n_quants = rowSums(!is.na(intensities.mtx) & intensities.mtx > 0.0),
                          is_full_quant = !is.na(rowSums(intensities.mtx > 0.0))) %>%
                  dplyr::rename_(.dots=intens_cols)

    message( 'Extracting pepmod states...' )
    pepmod_states.df <- dplyr::select(evidence.df, pepmod_state, pepmod_id, charge) %>% dplyr::distinct()

    # summarize intensities for pepmod_state X msrun pair (there could be multiple ones)
    message('Reshaping, summarizing & weighting pepmod_state intensities...')
    pms_intensities_long.df <- tidyr::gather(evidence.df, mstag, intensity, starts_with("intensity.")) %>%
      dplyr::mutate(mstag = factor(str_replace(mstag, "intensity\\.", ""),
                                   levels = levels(mschannels.df$mstag))) %>%
      dplyr::inner_join(mschannels.df %>% dplyr::select(msrun, mstag, mschannel)) %>%
      dplyr::group_by(pepmod_id, pepmod_state, msrun, mstag, mschannel) %>%
      dplyr::summarise(weight = weighted.mean(intensity, abs(1/mass_error_ppm)^0.5),
                       intensity = sum(intensity, na.rm=TRUE)) %>%
      dplyr::group_by(pepmod_id, mschannel) %>%
      dplyr::mutate(weight = pmax(weight/sum(weight), na_weight),
                    intensity = if_else(!is.na(intensity) & intensity>0, intensity, NA_real_)) %>%
      dplyr::ungroup()

    pms_full_intensities_long.df <- tidyr::expand(pms_intensities_long.df,
                                       mschannel, pepmod_state) %>%
      dplyr::left_join(pepmod_states.df) %>%
      dplyr::left_join(dplyr::select(mschannels.df, mschannel, msrun, mstag, one_of("msprotocol"))) %>%
      dplyr::left_join(pms_intensities_long.df) %>%
      dplyr::mutate(observed = !is.na(intensity) & intensity > 0,
                    weight = if_else(observed & is.finite(weight), weight, na_weight)) %>%
      dplyr::left_join(pepmod_states.df)
    if (!("msprotocol" %in% colnames(pms_full_intensities_long.df))) {
      pms_full_intensities_long.df$msprotocol <- "default"
    }

    message('Correcting & predicting missing intensities...')
    pms_full_intensities_long.df <- dplyr::group_by(
        dplyr::filter(pms_full_intensities_long.df, mstag != 'Sum') %>%
        dplyr::mutate(intensity_fixed = if_else(!is.na(intensity) & intensity > min_intensity, intensity, min_intensity)),
        pepmod_id, msprotocol) %>%
      dplyr::do({glm_corrected_intensities(., glm_max_factor = glm.max_factor, glm_reldelta_range = glm.reldelta_range)}) %>%
      dplyr::ungroup() %>%
      dplyr::select(-intensity_fixed) %>% # remove temporary non-NA column
      bind_rows(dplyr::filter(pms_full_intensities_long.df, mstag == 'Sum'))

    message('Summing intensities of different charges...')
    full_intensities_long.df <- dplyr::group_by(pms_full_intensities_long.df,pepmod_id, msrun) %>%
        dplyr::group_by(pepmod_id, msrun, mschannel, mstag) %>%
        dplyr::summarise(intensity = sum(intensity, na.rm=TRUE),
                         intensity_glm = sum(intensity_glm, na.rm=TRUE),
                         intensity_corr = sum(intensity_corr, na.rm=TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(intensity = if_else(intensity > 0, intensity, NA_real_),
                      intensity_glm = if_else(mstag == 'Sum', NA_real_, intensity_glm),
                      intensity_corr = if_else(intensity_corr > glm.min_intensity, intensity_corr, NA_real_),
                      intensity_glm_reldelta = (intensity_glm - intensity) / pmax(intensity_glm, intensity),
                      intensity_corr_reldelta = (intensity_corr - intensity) / pmax(intensity_corr, intensity))

    message( 'Extracting pepmod information...' )
    pepmodXmsrun_stats.df <- evidence.df %>% dplyr::group_by(pepmod_id, msrun) %>%
        dplyr::summarize(n_charges = n_distinct(charge),
                         n_evidences = n_distinct(evidence_id),
                         n_quants = sum(n_quants > 0),
                         n_full_quants = sum(is_full_quant)) %>%
        dplyr::ungroup()
    prediction.df <- dplyr::group_by(pms_full_intensities_long.df, pepmod_id) %>%
      # FIXME for different msprotocols GLM could be different
      dplyr::summarise(glm_rhs = glm_rhs[1]) %>% dplyr::ungroup() %>%
      dplyr::mutate(glm_rhs = factor(glm_rhs))
    pepmod_stats.df <- pepmodXmsrun_stats.df %>%
        dplyr::group_by(pepmod_id) %>%
        dplyr::summarise(n_msruns = n(),
                         n_quant_msruns = sum(n_quants > 0),
                         n_full_quant_msruns = sum(n_full_quants > 0),
                         n_max_charges=max(n_charges),
                         n_max_evidences=max(n_evidences)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(prediction.df)

    pepmods.df <- dplyr::select(evidence.df,
                                one_of(c("pepmod_id", "peptide_id", "protgroup_ids",
                                "protein_acs", "lead_protein_acs", "lead_razor_protein_ac",
                                "gene_names", "protein_names", "is_reverse", "is_contaminant",
                                "seq", "seq_len", "count_K", "count_R", "modifs", "n_miscleaved",
                                # Carbamidomethyl (C) site IDs', Oxidation (M) site IDs
                                "charge"))) %>%
                  dplyr::distinct() %>%
                  dplyr::group_by(pepmod_id) %>%
                  dplyr::mutate(charges = paste0(sort(charge), collapse=' ')) %>%
                  dplyr::filter(row_number() == 1L) %>%
                  dplyr::ungroup() %>%
                  dplyr::select(-charge) %>%
        left_join(pepmod_stats.df) %>%
        dplyr::mutate(is_shared_by_groups = grepl( ';', protgroup_ids, fixed = TRUE ),
                      is_shared = is_shared_by_groups,
                      is_shared_by_proteins = grepl( ';', protein_acs, fixed = TRUE ))

    message( 'Extracting peaks information...' )
    peak_columns <- c( 'pepmod_id', 'msrun', 'evidence_id', 'n_quants', 'is_full_quant',
    # Carbamidomethyl (C) Probabilities, Oxidation (M) Probabilities, Carbamidomethyl (C) Score Diffs, Oxidation (M) Score Diffs, Acetyl (Protein N-term), Carbamidomethyl (C), Oxidation (M),
      'ident_type', 'charge', 'pepmod_state', 'label_state', 'ms2_mz', 'mz', 'mass_da', 'resolution',
      'delta_mz_ppm', 'delta_mz_da', 'mass_error_ppm', 'Mass Error [mDa]', 'Mass Error [Da]',
      'Uncalibrated Mass Error [ppm]', 'Uncalibrated Mass Error [mDa]', 'Uncalibrated Mass Error [Da]',
      'Max intensity m/z 0', 'Max intensity m/z 1',
      'Retention time', 'rt_len', 'rt_avg', 'rt_start', 'rt_finish',
      'Retention time calibration', 'Match time difference', 'Match q-value', 'Match score',
      'Number of data points', 'Number of scans', 'Number of isotopic peaks', 'PIF', 'Fraction of total spectrum',
      'Base peak fraction', 'PEP', 'MS/MS Count', 'MS/MS Scan Number', 'Score', 'Delta score', 'Combinatorics',
      'MS/MS IDs', 'Best MS/MS', 'AIF MS/MS IDs' ) %>%
      .[ . %in% colnames(evidence.df) ]
    peaks.df <- evidence.df %>% dplyr::select_( .dots = backquote( peak_columns ) ) %>% dplyr::distinct()

    ratio_columns.df <- expand.grid( measure = 'ratio',
                                     mstag_nom = unique(dplyr::filter(intensity_columns.df, type == 'measured')$mstag),
                                     mstag_denom = unique(dplyr::filter(intensity_columns.df, type == 'measured')$mstag),
                                     type = c( '', 'normalized', 'shift' ) ) %>%
      dplyr::filter(("ratio" %in% import_data) & mstag_nom != mstag_denom) %>%
      mutate( old_name = gsub('\\s$', '', paste0('Ratio ', mstag_nom, '/', mstag_denom, ' ', type)),
              inverted_old_name = gsub('\\s$', '', paste0('Ratio ', mstag_denom, '/', mstag_nom, ' ', type)),
              suffix = paste0(type, if_else(type != "", '_', ""), mstag_nom, mstag_denom), 
              new_name = paste0( measure, '.', suffix ),
              exists = old_name %in% colnames(evidence.df),
              inverted_exists = inverted_old_name %in% colnames(evidence.df))

    # quant_mode defines which features are used for quantitation
    # label-free mode favours features that are observed in multiple runs
    # labeled mode favours features that have most quantitations per label

    if ("ratio" %in% import_data) {
      if ( correct_ratios ) {
        # channel/msrun ratios are more precisely measured than intensities even within the features of same pepmod
        # so calculate average ratios per pepmod and use these ratios to correct the absolute intensities:
        # the ratios of the corrected intensities should match
        message( "Correcting intensities by averaging ratios..." )
        if (quant_mode == "label-free") {
          if (!("condition" %in% colnames(mschannels.df))) {
            stop("No condition annotations, cannot average ratios, specify mschannel_annotate.f function")
          }
          rlm.control <- lmrob.control("KS2014")
          intensities.df <- dplyr::left_join(intensities.df, mschannels.df) %>%
              dplyr::left_join(dplyr::select(pepmods.df, pepmod_id, protgroup_ids)) %>%
              dplyr::group_by(condition, protgroup_ids)
          intensities.df <- dplyr::do(intensities.df, {glm_corrected_intensities(., rlm.control)})
          intensities.df <- intensities.df %>%
              dplyr::ungroup() %>% dplyr::select(-protgroup_ids, -condition)
        } else if (quant_mode == "labeled") {
          ref_ratio_cols.df <- bind_rows(
            dplyr::filter(ratio_columns.df, exists & mstag_denom == correct_by_ratio.ref_label) %>%
              dplyr::mutate(is_inverted = FALSE,
                            trf_old_name = old_name),
            dplyr::filter(ratio_columns.df, exists & mstag_nom == correct_by_ratio.ref_label) %>%
              dplyr::mutate(mstag_nom = mstag_denom,
                            mstag_denom = correct_by_ratio.ref_label,
                            is_inverted = TRUE,
                            trf_old_name = gsub( '\\s$', '', paste0( measure, ' ', mstag_nom, '/', mstag_denom, ' ', type )))) %>%
            dplyr::filter(type == "")
          inverted_cols <- backquote(ref_ratio_cols.df$old_name[ref_ratio_cols.df$is_inverted])
          names(inverted_cols) <- ref_ratio_cols.df$trf_old_name[ref_ratio_cols.df$is_inverted]
          pre_intensities.df$ref_intensity <- pre_intensities.df[[paste0("Intensity ", correct_by_ratio.ref_label)]]

          agg_ratios.df <- pre_intensities.df %>% dplyr::select(pepmod_id, msrun, ref_intensity, starts_with("Ratio"), mass_error_ppm) %>%
              dplyr::mutate(ratio_weight = abs(1/mass_error_ppm)/sum(abs(1/mass_error_ppm), na.rm=TRUE)) %>%
              dplyr::mutate_at(.cols = inverted_cols, funs(1/.)) %>%
              dplyr::summarise_each_(funs(if_else(all(is.na(ratio_weight)), NA_real_,
                                               weighted.mean(.[!is.na(ratio_weight)],
                                                              ratio_weight[!is.na(ratio_weight)], na.rm=TRUE))),
                                     backquote(ref_ratio_cols.df$trf_old_name))
          ilabels_ordered <- c(as.character(ref_ratio_cols.df$mstag_nom), correct_by_ratio.ref_label)
          intens_mtx <- as.matrix(intensities.df[paste0("intensity.", ilabels_ordered)])
          ratios_mtx <- cbind(as.matrix(agg_ratios.df[ref_ratio_cols.df$trf_old_name]),
                              if_else(is.na(intens_mtx[[ncol(intens_mtx)]]), 0, 1))
          w_avg_intens <- rowSums(intens_mtx * ratios_mtx, na.rm = TRUE)
          ratio_sqr_sum <- rowSums(ratios_mtx * ratios_mtx, na.rm = TRUE)
          # corrected reference intensity that would minimize
          # the total sqr deviation of corrected intensities vs original intensities
          ref_intensity_corr <- w_avg_intens/ratio_sqr_sum
          intens_corr_mtx <- ratios_mtx * replicate(ncol(ratios_mtx), ref_intensity_corr)
          intensities.df[paste0("intensity_corr.", ilabels_ordered)] <- intens_corr_mtx
          intensities.df$intensity_corr.Sum <- rowSums(intens_corr_mtx, na.rm = TRUE)
          # restore intensities that are NA because no ratio available
          for (lbl in ilabels) {
            col_corr <- paste0("intensity_corr.", lbl)
            col_orig <- paste0("intensity.", lbl)
            intensities.df[col_corr] <- if_else(is.na(intensities.df[[col_corr]]) | intensities.df[[col_corr]]==0.0,
                                                intensities.df[[col_orig]],
                                                intensities.df[[col_corr]])
          }
        }
    }
    }
    message( 'Converting intensities to ', match.arg(layout), ' row format...' )
    if ( match.arg(layout) == 'pepmod_mschannel' ) {
        message( 'Converting intensities to pepmod_mschannel row format...' )
        #print( str( intensities.df ) )
        intensities.df <- full_intensities_long.df %>% dplyr::filter(observed) %>%
          dplyr::mutate( mstag = factor( mstag, levels = as.character(intensity_columns.df$mstag) ) )
    } else if (match.arg(layout) == 'pepmod_msrun') {
        fullest.df <- tidyr::gather(full_intensities_long.df, mod, val, starts_with("intensity")) %>%
          dplyr::mutate(mod = str_replace(mod, "intensity", ""),
                        intensity = interaction(mod, mstag, lex.order = TRUE, sep='.')) %>%
          dplyr::filter(mstag != 'Sum' | mod == "") %>%
          dplyr::select(-mod, -mschannel, -mstag)
        intensities.df <- tidyr::spread(fullest.df, intensity, val, sep = "")
    }
    ident_types.df <- dplyr::distinct(dplyr::select(peaks.df, pepmod_id, msrun, ident_type)) %>%
      dplyr::mutate(has_ident = TRUE) %>%
      tidyr::spread(ident_type, has_ident, sep=".")
    intensities.df <- dplyr::left_join(intensities.df, ident_types.df) %>%
      mutate_at(.cols=vars(starts_with("ident_type")), funs(!is.na(.)))

    ratio_columns_sel.df <- dplyr::filter(ratio_columns.df,
                                          mstag_nom < mstag_denom & mstag_denom != 'Sum')
    if (any(ratio_columns_sel.df$exists)) {
      message('Extracting ratio data...')
      if (!all(ratio_columns_sel.df$exists)) {
        missing_cols <- dplyr::filter(ratio_columns_sel.df, !exists) %>% .$old_name
        warning(nrow(missing_cols), " ratio column(s) are missing: ", paste(missing_cols, collapse=' '))
      }
      sel_cols <- backquote(ratio_columns_sel.df$old_name)
      names(sel_cols) <- ratio_columns_sel.df$new_name
      ratios.df <- intensities.df %>%
                  dplyr::summarise_at(sel_cols, funs(mean(., na.rm=TRUE))) %>%
                  dplyr::ungroup()
      ratios.df[,dplyr::filter(ratio_columns_sel.df,!exists)$new_name] <- NA_real_ # add missing column names to make it evidence.df hapeable
      # fix NA
      ratios.df <- ratios.df %>% mutate_at(ratio_columns_sel.df$new_name,
                                           funs(if_else(!is.na(.) & .==0, NA_real_, .)))
      if (match.arg(layout) == 'long') {
          message( 'Converting ratios to long format...' )
          ratios.df <- reshape( ratios.df, direction = 'long',
                  idvar = id_cols,
                  #varying = ratio_columns_sel.df$new_name,
                  v.names = 'ratio',
                  timevar = 'ratio_type',
                  times = ratio_columns_sel.df$suffix,
                  sep = '.') %>%
            dplyr::inner_join(dplyr::select(ratio_columns_sel.df, suffix, type, mstag_nom, mstag_denom),
                              by = c("ratio_type" = "suffix"))
      }
    } else {
      message( 'No channel ratios')
      ratios.df <- NULL
    }
    protgroup_ids <- pepmods.df$protgroup_ids %>% unique()
    protgroups2protgroup.list <- strsplit(protgroup_ids, ';', fixed=TRUE)
    protgroups2protgroup.df <- data.frame(protgroup_ids = rep.int(protgroup_ids, sapply(protgroups2protgroup.list, length)),
                                          protgroup_id = as.integer(unlist(protgroups2protgroup.list)),
                                          stringsAsFactors = FALSE )
    return ( list( pepmods = pepmods.df,
                   protgroups2protgroup = protgroups2protgroup.df,
                   raw_files = msruns.df,
                   mschannels = mschannels.df,
                   peaks = peaks.df,
                   intensities = intensities.df,
                   ratios = ratios.df ) )
}

read.MaxQuant.Evidence <- function( folder_path, file_name = 'evidence.txt', layout = c( "wide", "long" ), nrows = -1,
                                    nrows = Inf, guess_max = min(10000L, nrows),
                                    min_pepmod_state_freq = 0.9, min_essential_freq = 0.0,
                                    correct_ratios = TRUE, correct_by_ratio.ref_label = NA,
                                    mode = c("labeled", "label-free"), mschannel_annotate.f = NULL )
{
  process.MaxQuant.Evidence(read.MaxQuant.Evidence_internal(
                            folder_path = folder_path, file_name = file_name, nrows = nrows, guess_max=guess_max),
                            min_pepmod_state_freq = min_pepmod_state_freq,
                            min_essential_freq = min_essential_freq,
                            correct_ratios = correct_ratios,
                            correct_by_ratio.ref_label = correct_by_ratio.ref_label,
                            mode = mode, mschannel_annotate.f = mschannel_annotate.f)
}

msrun_code_parser.f = function( msruns, chunk_names = c( 'dataset', 'batch', 'frac_protocol', 'fraction', 'tech_replicate' ) ) {
      msrun_chunks <- strsplit( as.character( msruns ), '_' )
      n_missing_chunks <- length( chunk_names ) - sapply( msrun_chunks, length )
      msrun_chunks <- lapply( seq_along(msrun_chunks), function(i) c(msrun_chunks[[i]],rep.int(NA,n_missing_chunks[[i]])) )
      res <- cbind( msrun = msruns, as.data.frame( do.call( rbind, msrun_chunks ), stringsAsFactors = FALSE ) )
      colnames(res) <- c( 'msrun', chunk_names )
      canonical_order <- c('msrun','dataset','experiment','batch','frac_protocol','replicate','fraction','tech_replicate')
      res[,setdiff(canonical_order,colnames(res))] <- NA
      res[,canonical_order]
}

expand_protgroup_acs <- function(protgroups, acs_col, ac_col = str_replace(acs_col, "_acs", "_ac") ) {
  acs <- str_split(protgroups[[acs_col]], ";", simplify = FALSE)
  res <- data.frame(protgroup_id = rep.int(protgroups$protgroup_id, sapply(acs, length)),
                    stringsAsFactors = FALSE)
  res[[ac_col]] <- unlist(acs)
  return (res)
}

match_protgroup_by_acs <- function(pgs1, pgs2, acs_col, suffix=c(".x", ".y")) {
  ac_col = str_replace(acs_col, "_acs", "_ac")
  expd_pgs1 = expand_protgroup_acs(pgs1, acs_col, ac_col)
  expd_pgs2 = expand_protgroup_acs(pgs2, acs_col, ac_col)
  dplyr::full_join(expd_pgs1, expd_pgs2, by=ac_col, suffix=suffix)
}

split_protgroups <- function(protgroups.df) {
    majority_acs.df <- expand_protgroup_acs(protgroups.df, acs_col = 'majority_protein_acs') %>%
        dplyr::mutate(is_majority = TRUE)
    expand_protgroup_acs(protgroups.df, acs_col = 'protein_acs') %>%
    mutate( is_contaminant = str_detect(protein_ac, "^CON_"),
            is_reverse = str_detect(protein_ac, "^REV_") ) %>%
    dplyr::left_join(majority_acs.df, by=c("protgroup_id"="protgroup_id", "protein_ac" = "majority_protein_ac")) %>%
      dplyr::mutate(is_majority = if_else(is.na(is_majority), FALSE, TRUE))
}
