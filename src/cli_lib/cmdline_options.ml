(* Common arguments shared between programs. *)

module Ns = String
open Prohlatype
open Cmdliner

let (//) = Filename.concat

let repo = "prohlatype"
let bugs =
  sprintf "Browse and report new issues at <https://github.com/hammerlab/%s>"
        repo

let author = "Leonid Rozenberg <leonidr@gmail.com>"

let version = Version.version

(*** Basic command line parsers and printers. ***)
let positive_int_parser w s =
  try
    let d = Scanf.sscanf s "%d" (fun x -> x) in
    if d <= 0 then
      Error (`Msg (s ^ " is not a positive integer"))
    else
      Ok (w d)
  with Scanf.Scan_failure msg ->
    Error (`Msg msg)

let int_fprinter frmt =
  Format.fprintf frmt "%d"

let positive_int =
  Arg.conv ~docv:"POSITIVE INTEGER"
    ((positive_int_parser (fun n -> n)), int_fprinter)

let positive_float_parser s =
    try
      let f = Scanf.sscanf s "%f" (fun x -> x) in
      if f <= 0. then
        Error (`Msg (s ^ " is not a positive float"))
      else
        Ok f
    with Scanf.Scan_failure msg ->
      Error (`Msg msg)

let positive_float =
  Arg.conv ~docv:"POSITIVE FLOAT"
    (positive_float_parser, fun frmt -> Format.fprintf frmt "%f")

let greater_than_one =
  let docv = "natural number greater than 1" in
  Arg.conv ~docv:(Ns.uppercase_ascii docv)
    ((fun s ->
      try
        let d = Scanf.sscanf s "%d" (fun x -> x) in
        if d <= 1 then
          Error (`Msg (s ^ " is not a " ^ docv))
        else
          Ok d
      with Scanf.Scan_failure msg ->
        Error (`Msg msg)), int_fprinter)

let non_negative_int_parser =
  fun s ->
    try
      let d = Scanf.sscanf s "%d" (fun x -> x) in
      if d < 0 then
        Error (`Msg (s ^ " is negative"))
      else
        Ok d
    with Scanf.Scan_failure msg ->
      Error (`Msg msg)

let non_negative_int =
  Arg.conv ~docv:"NON-NEGATIVE INTEGER"
    (non_negative_int_parser, int_fprinter)

(*** Graph source arguments. ***)

(* Either a well formed path to an alignment file or a well formed path prefix
 * to valid cDNA (eg. A_nuc.txt) and gDNA (eg. A_gen.txt) is necessary for
 * typing, graphing, making sequences ... etc.
 *)
let error_msg fmt =
  ksprintf (fun s -> Error (`Msg s)) fmt

let recognized_locus path =
  let base = Filename.basename path in
  let wext =
    try Filename.chop_extension base
    with (Invalid_argument _) -> base
  in
  match String.split ~on:(`Character '_') wext with
  | [pl]
  | [pl; "gen"]
  | [pl; "nuc"] ->
    begin match Nomenclature.parse_locus pl with
    | Ok o      -> Ok o
    | Error e   -> error_msg "%s" e
    end
  | _           -> error_msg "Unrecognized file or prefix: %s" path

(*** Allele selector arguments. ***)
let regex_command_line_arg = "allele-regex"
let allele_command_line_arg = "specific-allele"
let number_command_line_arg = "number-of-alleles"
let do_not_ignore_suffixed_alleles_flag = "do-not-ignore-suffixed-alleles"

let allele_selector_paragraph reason = sprintf
  "Various command line options allow the user to pare down, \
  and construct the set of alleles that are included, via \
  \"allele-selectors\". This way, %s two alleles (eg. A*01:01:01:01 \
  vs A*02:01:01:01) or two branches of the allele tree (eg. A*01 vs A*03). \
  Initially, the set of alleles consists of all the ones found in the \
  alignment file except for those with a suffix (see --%s). Afterwards if we \
  find one allele-selector, the set is reset to empty and alleles are added \
  by applying each allele selectors."
  reason
  do_not_ignore_suffixed_alleles_flag

let regex_arg =
  let docv = "REGEX" in
  let doc  = sprintf
    "An allele-selector to specify alleles to add to the allele-set via a \
     regex. This option is similar to --%s, but lets the user specify a wider \
     range (specifically those alleles matching a POSIX regex) of alleles to \
     include.  The '*' used in HLA allele names must be properly escaped from \
     the command line: i.e. --%s \"A\\\\*02:03\". \ One can repeat this \
     argument to specify different groups of alleles. "
    allele_command_line_arg
    regex_command_line_arg
  in
  let open Arg in
  let parser_ = parser_of_kind_of_string ~kind:docv
    (fun s -> Some (Alleles.Selectors.Regex s))
  in
  value & opt_all (conv ~docv (parser_, Alleles.Selectors.pp)) []
  & info ~doc ~docv [regex_command_line_arg]

let specific_arg =
  let docv = "STRING" in
  let doc  = sprintf
    "Specify specific alleles to add to the allele-set. \
     One can repeat this argument to specify multiple alleles."
  in
  let open Arg in
  let parser_ =
    parser_of_kind_of_string ~kind:docv
      (fun s -> Some (Alleles.Selectors.Specific s))
  in
  value & opt_all (conv ~docv (parser_, Alleles.Selectors.pp)) []
        & info ~doc ~docv [allele_command_line_arg]

let number_alleles_arg =
  let docv = "POSITIVE INTEGER" in
  let doc  = "Number of alleles to add to the allele-set. \
              This allele-selectors can be to resize the number of alleles \
              considered from a given gene." in
  let open Arg in
  let parser_ = (positive_int_parser (fun d -> Alleles.Selectors.Number d)) in
  let nconv = conv ~docv (parser_, Alleles.Selectors.pp) in
  (value & opt (some nconv) None & info ~doc ~docv [number_command_line_arg])

let do_not_ignore_suffixed_alleles_flag =
  let doc =
    "IMGT publishes alleles with suffixes ('N', 'L', 'S', 'C', 'A' or 'Q') to \
     indicate special conditions such as null alleles ('N') or questionable \
     expression ('Q'). Often times, these alleles can complicate prohlatype as \
     reads sometimes align surprisingly well to them (eg. HLA-A*11:50Q) as \
     opposed to the true allele. We want to assign these a small prior \
     probability, smaller than the other alleles. For the purposes of having a \
     clearer typing algorithm we exclude all such alleles from the starting \
     allele-set by default. For rare cases or diagnostic purposes one can pass \
     this flag to bring them back into the starting allele set."
  in
  Arg.(value & flag & info ~doc [do_not_ignore_suffixed_alleles_flag])

(*** Other args. ***)
let no_cache_flag item =
  let doc =
    sprintf "Do not use a disk cache (in %s subdirectory of the current \
             directory) to search for previously constructed %s. \
             This flag will also turn off the default behavior of saving \
             this %s there."
      Cache.dir item item
  in
  Arg.(value & flag & info ~doc ["no-cache"])

let do_not_join_same_sequence_paths_flag =
  let doc = "Do not join same sequence paths; remove overlapping nodes at the \
             same position, in the string graph."
  in
  Arg.(value & flag & info ~doc ["do-not-join-same-sequence-paths"])

let option_to_list o =
  Option.value_map o ~default:[] ~f:(fun s -> [s])

let aggregate_selectors
    ~number_alleles
    ~regex_list
    ~specific_list
    ~do_not_ignore_suffixed_alleles =
    let acc =
      regex_list
      @ specific_list
      @ (option_to_list number_alleles)
    in
    if do_not_ignore_suffixed_alleles then
        Alleles.Selectors.DoNotIgnoreSuffixed :: acc
    else
      acc

let verbose_flag =
  let doc = "Print progress messages to stdout." in
  Arg.(value & flag & info ~doc ["v"; "verbose"])

let kmer_size_arg =
  let default = 10 in
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Number of consecutive nucleotides to use consider in K-mer
              index construction. Defaults to %d." default
  in
  Arg.(value & opt positive_int default & info ~doc ~docv ["k"; "kmer-size"])

let fastq_file_arg =
  let docv = "FASTQ-FILE" in
  let doc = "Fastq formatted DNA reads file. Only one file per sample. \
             List paired end reads as 2 sequential files." in
  Arg.( non_empty
      & pos_right 0 file []
      & info ~doc ~docv [])

let num_reads_arg =
  let docv = "POSITIVE INTEGER" in
  let doc =
    "Number of reads to use from the front of the FASTA file. \
    This option is intended to stop execution after a given number of samples."
  in
  Arg.(value & opt (some positive_int) None & info ~doc ~docv ["reads"])

let trie_argument = "trie"
let weighted_per_segment_argument = "weighted-segment"
let reference_distance_argument = "reference-distance"

let optional_distance_flag, defaulting_distance_flag =
  let open Distances in
  let opts =
    [ Trie
      , "trie based off of allele names"
      , false
      , trie_argument
    ; WeightedPerSegment
      , "smallest shared weighted per segment distance"
      , true
      , weighted_per_segment_argument
    ; Reference
      , "consider the reference the closest for all alleles"
      , false
      , reference_distance_argument
    ]
  in
  let r default =
    fun s ->
    sprintf "How to compute the distance between alleles: %s.%s"
      s (if default then "default" else "")
  in
  let o = sprintf "How to compute the distance between alleles: %s. \
             Specifying a distance argument will turn on imputation for
             alleles over their missing segments."
  in
  let open Arg in
  value & vflag None
    (List.map opts ~f:(fun (opt,m,_,c) -> (Some opt, info ~doc:(o m) [c])))
  , value & vflag WeightedPerSegment
      (List.map opts ~f:(fun (opt,m,d,c) -> (opt, info ~doc:(r d m) [c])))

let print_top_flag =
  let doc = "Print only the specified number (positive integer) of alleles" in
  Arg.(value & opt (some int) None & info ~doc ["print-top"])

let specific_read_args =
  let docv = "STRING" in
  let doc  = "Read name string (to be found in FASTQ) to type. Add multiple \
              to create a custom set of reads." in
  Arg.(value
      & opt_all string []
      & info ~doc ~docv ["sr"; "specific-read"])

let reduce_resolution_arg =
  let doc  = "Reduce the resolution of the PDF, to a lower number of \
              \"digits\". The general HLA allele nomenclature \
              has 4 levels of specificity depending on the number of colons \
              in the name. For example A*01:01:01:01 has 4 and A*01:95 \
              has 2. This argument specifies the number (1,2,or 3) of \
              digits to reduce results to. For example, specifying 2 will \
              choose the best (depending on metric) allele out of \
              A*02:01:01:01, A*02:01:01:03, ...  A*02:01:02, ... \
              A*02:01:100 for A*02:01. The resulting set will have at most \
              the specified number of digit groups, but may have less."
  in
  let one_to_three_parser s =
    match positive_int_parser (fun x -> x) s with
    | Ok x when x = 1 || x = 2 || x = 3 -> Ok x
    | Ok x                              -> Error (`Msg (sprintf  "not 1 to 3: %d" x))
    | Error e                           -> Error e
  in
  let open Arg in
  let one_to_three = conv ~docv:"ONE,TWO or THREE"
    (one_to_three_parser , (fun frmt -> Format.fprintf frmt "%d"))
  in
  value & opt (some one_to_three) None & info ~doc ["reduce-resolution"]

let probability_parser s =
  try
    let f = Scanf.sscanf s "%f" (fun x -> x) in
    if f < 0.0 then
      Error (`Msg (s ^ " is less than zero, not a valid probability."))
    else if f = 0.0 then begin
      eprintf "zero probability specified; things might be strange!\n";
      Ok f
    end else if f = 1.0 then begin
      eprintf "one probability specified; things might be strange!\n";
      Ok f
    end else if f > 1.0 then
      Error (`Msg (s ^ " is greater than one, not a valid probability."))
    else
      Ok f
  with Scanf.Scan_failure msg ->
    Error (`Msg msg)

let probability_arg =
  Arg.conv ~docv:"PROBABILITY"
    (probability_parser, (fun fmt -> Format.fprintf fmt "%f"))

let insert_probability_arg =
  let docv = "PROBABILITY" in
  let doc  =
    sprintf "Specify a value between 0 and 1 to represent the insert \
              emission probability. The default value is: %f"
      Phmm.default_insert_probability
  in
  Arg.(value
      & opt probability_arg Phmm.default_insert_probability
      & info ~doc ~docv ["insert-emission-probability"])

let max_number_mismatches_arg =
  let docv = "POSITIVE INTEGER" in
  let doc = "Setup a filter on the reads to cancel evaluation once we've \
             seen this many mismatches.
             This will restrict which alleles are used for inference, \
             but in aggregate take less time." in
  Arg.(value & opt (some positive_int) None & info ~doc ~docv
        ["max-mismatches"])

let read_size_override_argument = "read-size"

let read_size_override_arg =
  let docv = "POSITIVE INTEGER" in
  let doc = "Override the number of bases to calculate the likelihood over, \
              instead of using the number of bases in the FASTQ." in
  Arg.(value
      & opt (some positive_int) None
      & info ~doc ~docv [read_size_override_argument])

let allele_depth_info_argument = "allele-depth-info"
let allele_depth_info_default = 5

let allele_depth_info_arg =
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Specify a positive integer to indicate the number of best alleles \
             and positions to report. Defaults to %d. Be cautious about \
             increasing this number as keeping track of this data can slow \
             down the final analysis."
      allele_depth_info_default
  in
  Arg.(value & opt positive_int allele_depth_info_default
             & info ~doc ~docv [allele_depth_info_argument])

let output_format_flag =
  Arg.(value & vflag `TabSeparated
    [ `TabSeparated
      , info ~doc:("Output results in a tab separated format. Default.") ["tab"]
    ; `Json
      , info ~doc:("Output results in JSON. Not the default.") ["json"]
    ])

(* Band arguments
let not_band_flag =
  let doc = "Calculate the full forward pass matrix." in
  Arg.(value & flag & info ~doc ["do-not-band"])
*)

let band_warmup_argument = "band-warmup"
let band_number_argument = "band-number"
let band_radius_argument = "band-radius"

let band_warmup_arg =
  let default = ParPHMM.(band_default.warmup) in
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "At which row (bases of the read) in the forward pass to start \
      computing bands instead of the full pass. Defaults to: %d.\n Passing \
      this (or any of the other two band settings %s, %s) argument will \
      trigger the banding logic."
      default (Arg.doc_quote band_number_argument) band_radius_argument
  in
  Arg.(value
        & opt (some positive_int) None
        & info ~doc ~docv [band_warmup_argument])

let number_bands_arg =
  let default = ParPHMM.(band_default.number) in
  let docv = "POSITIVE INTEGER" in
  let doc  =
    sprintf "Number of bands to calculate. Defaults to %d.\n When we perform \
      banded passes there is a chance that after the \"warmup\" period we \
      haven't correctly located where the read aligns to the reference; this \
      may easily occur if there is an error in first part of the read. To \
      compensate we can calculate more than one band and average over the \
      results as one of them will (probably) capture most of the probability \
      mass. Passing this (or any of the other two band settings %s, %s) \
      argument will trigger the banding logic."
      default band_warmup_argument band_radius_argument
  in
  Arg.(value
        & opt (some positive_int) None
        & info ~doc ~docv [band_number_argument])

let band_radius_arg =
  let default = ParPHMM.(band_default.radius) in
  let docv = "INTEGER greater than 1" in
  let doc  =
    sprintf "Radius of bands to calculate. Defaults to %d. Must be greater \
      than 1.\n For a given read position, bands have this specified radius \
      around the most likely match likelihood; the calculation logic tries to \
      calculate this many cells before and after in the reference position. \
      Passing this (or any of the other two band settings %s, %s) argument \
      will trigger the banding logic."
      default band_warmup_argument band_number_argument
  in
  Arg.(value
        & opt (some greater_than_one) None
        & info ~doc ~docv [band_radius_argument ])

let forward_pass_accuracy_arg =
  let docv = "POSITIVE DOUBLE" in
  let doc  = sprintf
    "Tweak the accuracy (and consequently speed) of the forward pass \
    (larger -> less accurate but faster). The default value is %f, \
    which is sufficient for most Phred error profiles.
    The implementation uses 63 bit floats which is more than sufficient for \
    this calculation."
    !Probability.dx
  in
  Arg.(value & opt (some positive_float) None
             & info ~doc ~docv ["numerical-accuracy"])

let do_not_past_threshold_filter_flag =
  let doc = "Do not use previously calculated likelihoods (from other loci \
             or the reverse complement alignment) as a threshold filter to \
             short-circuit evaluations."
  in
  Arg.(value & flag & info ~doc ["do-not-past-threshold-filter"])

let do_not_finish_singles_flag =
  let doc = "Paired FASTQ files that may be created by downstream analysis may \
             not guarantee that all of the reads are paired. When in paired \
             mode (2 FASTQ files are passed), by default we'll continue \
             evaluating the unpaired reads as if they were single. By passing \
             this flag the program will skip them entirely and use only the \
             paired reads."
  in
  Arg.(value & flag & info ~doc ["do-not-finish-singles"])

let likelihood_report_size_argument = "likelihood-report-size"
let likelihood_report_size_arg =
  let open ParPHMM_drivers.Output in
  let docv = "POSITIVE INTEGER" in
  let doc = "Override the default number of individual allele likelihoods to \
    report. In general, this is $(b,not) the desired outcome; we should \
    consider the likelihoods of diploids (referred to as zygosity here). \
    None-the-less this can be useful for relative allele diagnosis. By \
    default, a value will be reported for all alleles."
  in
  Arg.(value & opt (some non_negative_int) None
             & info ~doc ~docv [likelihood_report_size_argument])

let non_negative_zygosity =
  let open ParPHMM_drivers.Zygosity_best in
  let non_negative_int_parser =
    fun s ->
      try
        let d = Scanf.sscanf s "%d" (fun x -> x) in
        if d < 0 then
          Error (`Msg (s ^ " is negative"))
        else if d = 0 then
          Ok NoSpec
        else
          Ok (Spec d)
      with Scanf.Scan_failure msg ->
        Error (`Msg msg)
  in
  Arg.conv ~docv:"NON-NEGATIVE INTEGER"
    (non_negative_int_parser
    , fun frmt zas ->
        match zas with
        | NoSpec    -> Format.fprintf frmt "%d" 0
        | Spec d    -> Format.fprintf frmt "%d" d
        | NonZero v -> Format.fprintf frmt "P > %f" v)

let zygosity_report_size_argument = "zygosity-report-size"
let zygosity_report_size_arg =
  let open ParPHMM_drivers in
  let docv = "POSITIVE INTEGER" in
  let doc = sprintf
    "Override the default number of allelic pairs reported as part of \
    the zygosity portion. By defaults only the pairs that have non-zero \
    (> %f) probability will be reported. Set this value to zero to report \
    all (could be > 10e6) values."
    Zygosity_best.default_non_zero
  in
  Arg.(value & opt non_negative_zygosity
                      Zygosity_best.(NonZero default_non_zero)
             & info ~doc ~docv [zygosity_report_size_argument])

let zygosity_non_zero_value_arg =
  let open ParPHMM_drivers in
  let docv = "POSITIVE FLOAT" in
  let doc = sprintf
    "Override the default lowerbound of non-zero zygosities. The \
    default (%f) probability might be too high for some scenarios, such \
    as if there are too few reads. This argument overrides --%s."
    Zygosity_best.default_non_zero zygosity_report_size_argument
  in
  Arg.(value & opt (some positive_float) None
             & info ~doc ~docv ["zygosity-non-zero"])

let to_num_zygosities ~zygosity_non_zero_value ~zygosity_report_size =
  Option.value_map zygosity_non_zero_value
    ~f:(fun v -> ParPHMM_drivers.Zygosity_best.NonZero v)
    ~default:zygosity_report_size

let per_reads_report_size_argument = "per-read-report-size"
let per_reads_report_size_arg =
  let docv = "POSITIVE INTEGER" in
  let doc = sprintf
    "Override the number of (per)-read information to report. By \
    default all of the per-read information is reported."
  in
  Arg.(value & opt (some non_negative_int) None
             & info ~doc ~docv [per_reads_report_size_argument])

let number_of_processors_arg =
  let docv = "POSITIVE INTEGER" in
  let doc = sprintf
    "Parallelize read processing across this number of processesors. \
     By default the application is not parallelized as the user would \
     most likely achieve greater efficiency by naively parallelizing \
     across samples as opposed to parallelizing across reads. \
     Furthermore the user MUST specify --%s in this mode."
    read_size_override_argument
  in
  Arg.(value & opt (some positive_int) None
             & info ~doc ~docv ["number-of-processors"])

let split_argument = "split"
let split_arg =
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Split the processing of a read into this many smaller forward \
      passes, that are then aggregated. This should decrease the run time cost \
      of a full pass as the match probabilities, across alleles, are less \
      distributed. This value $(b,must) divide the read\ size evenly. \
      Be careful about splitting too much (setting this value too high) as \
      longer reads align more accurately than shorter reads (specifically, \
      10bp is too short while 25bp rarely gets filtered). If the resulting \
      output results in too many filtered reads, consider using a smaller \
      value."
  in
  Arg.(value & opt (some positive_int) None
             & info ~doc ~docv [ split_argument ])

let not_prealigned_flag =
  let doc = "There are 2 PHMM implementations that differ in what kind of \
             transitions they allow. The default assumes that reads will fit \
             completely within the gene reference region (this allows us to \
             make small optimizations). Pass this flag if the reads are not \
             guaranteed to have this property and may overlap partly with the \
             start or end of the region."
  in
  Arg.(value & flag & info ~doc ["not-prealigned"])

let output_arg =
  let docv = "FILENAME PREFIX" in
  let doc = sprintf
    "Change output from stdout to files that are prefixed with the \
    specified argument. A log file will be generated with the \".log\" \
    suffix. Data output will be written to a file with either the \".tsv\" \
    or \".json\" suffix (see --json)."
  in
  Arg.(value & opt (some string) None
             & info ~doc ~docv ["o"; "output"])

let setup_oc output format_ =
  match output with
  | None   -> stdout, stdout
  | Some f ->
      let suffix =
        match format_ with
        | `TabSeparated -> "tsv"
        | `Json -> "json"
      in
      register_oc (sprintf "%s.log" f)
      , register_oc (sprintf "%s.%s" f suffix)


let genetic_argument = "genetic"
let nuclear_argument = "nuclear"
let merged_argument = "merged"

let class1_directory_argument = "class1"
let full_class1_directory_argument = "full-class1"

let gen_nuc_merged_flag =
  let gen_doc =
    "When a directory is specified, select the genetic files, eg. A_gen.txt. \
     This is the default choice."
  in
  let nuc_doc =
    "When a directory is specified select the nuclear files, eg. A_nuc.txt."
  in
  let mgd_doc =
    sprintf
      "When a directory is specified merge the genetic (eg. A_gen.txt) with \
       nuclear files (A_nuc.txt). "
  in
  Arg.(value & vflag `Genetic
         [ `Genetic, info ~doc:gen_doc [genetic_argument]
         ; `Nuclear, info ~doc:nuc_doc [nuclear_argument]
         ; `Merged,  info ~doc:mgd_doc [merged_argument]
         ])

let anded lst =
  let n = List.length lst in
  let before, after = List.split_n lst (n - 1) in
  let b = String.concat ~sep:", " before in
  let a = List.hd_exn after in
  sprintf "%s and %s" b a


let to_files gen_nuc_mgd d lg =
  let suffix =
    match gen_nuc_mgd with
    | `Genetic  -> "_gen.txt"
    | `Nuclear  -> "_nuc.txt"
    | `Merged   -> ""
  in
  Nomenclature.locus_group_to_loci lg
  |> List.map ~f:(fun l -> d // (Nomenclature.show_locus l ^ suffix))

(* These are convenience selectors that choose more than one locus at a time. *)
let class_selectors class1 full_class1 gen_nuc_mgd alignment_files merge_files =
  let files =
    match full_class1, class1 with
    | Some d,          _      -> to_files gen_nuc_mgd d Nomenclature.FullClassI
    | _     ,          Some d -> to_files gen_nuc_mgd d Nomenclature.ClassI
    | None,            None   -> []
  in
  if gen_nuc_mgd = `Merged then
    alignment_files, files @ merge_files
  else
    files @ alignment_files, merge_files

let class1_loci_str =
  let open Nomenclature in
  anded (List.map ~f:show_locus (locus_group_to_loci ClassI))

let full_class1_loci_str =
  let open Nomenclature in
  anded (List.map ~f:show_locus (locus_group_to_loci FullClassI))

let only_class1_flag = "only-common-class1"

let only_class1_directory_flag =
  let open Nomenclature in
  let doc = sprintf
    "Consider only the common class I HLA loci (%s) in the passed directory."
      class1_loci_str
  in
  Arg.(value & flag & info ~doc [only_class1_flag])

(* fpd - file, prefix or directory *)
let fpd_parser is_file_k path =
  if Sys.file_exists path then begin
    is_file_k path
  end else
    recognized_locus path >>= function
      | Nomenclature.P  ->
        let g = path ^ "_gen.txt" in
        if Sys.file_exists g then begin
          eprintf "As of 2017-12-15 P_nuc.txt doesn't exist. We'll use just \
            the genetic file.\n";
          Ok (`Prefix path)
        end else
          error_msg "Genetic alignment file doesn't exist: %s" g
      | l               ->
        if not (Alter_MSA.supported_loci l) then
          error_msg "Merging of gene %s not supported."
            (Nomenclature.show_locus l)
        else
          let n = path ^ "_nuc.txt" in
          let g = path ^ "_gen.txt" in
          if not (Sys.file_exists n) then begin
            error_msg "Nuclear alignment file doesn't exist: %s" n
          end else if not (Sys.file_exists g) then
            error_msg "Genetic alignment file doesn't exist: %s" g
          else
            Ok (`Prefix path)


let fp_printer frmt s = match s with
  | `File s      -> Format.fprintf frmt "File: %s" s
  | `Prefix s    -> Format.fprintf frmt "Prefix: %s" s

let fp_constructor =
  let is_locus path =
    recognized_locus path >>= fun l -> Ok (`File path)
  in
  Arg.conv ~docv:"File, prefix or directory"
    (fpd_parser is_locus, fp_printer)

let file_prefix_doc ~action ~dir_str = sprintf
  "Determine what IMGT/HLA alignment data source to %s. \
  If a file is specifed (eg. path-to-IMGT/alignments/A_gen.txt) then \
  only the alleles in that file are used. \
  %s
  If a path to a $(b,file prefix) is specified \
  (eg. path-to-IMGT/alignments/A) then a \"A_gen.txt\" and \"A_nuc.txt\" \
  file must also exist in the same directory and the sequence \
  information of the alleles in those files are merged. \
  See the imputation and distance (eg. --%s) logic for further \
  information."
  action
  dir_str
  weighted_per_segment_argument

let file_or_prefix_arg ~action =
  let docv = "FILE|PATH-to-PREFIX" in
  let doc = file_prefix_doc ~action ~dir_str:"" in
  Arg.(required
      & pos 0 (some fp_constructor) None
      & info ~doc ~docv [])

let file_or_prefix_to_allele_input ?distance ~selectors =
  let open Alleles.Input in function
  | `File path   -> alignment ?distance ~selectors path
  | `Prefix p    -> merge ?distance ~selectors p

let fpd_printer frmt s = match s with
  | `File s      -> Format.fprintf frmt "File: %s" s
  | `Prefix s    -> Format.fprintf frmt "Prefix: %s" s
  | `Directory s -> Format.fprintf frmt "Directory: %s" s

let fpd_constructor =
  let check_dir path =
    if Sys.is_directory path then
      Ok (`Directory path)
    else
      recognized_locus path >>= fun l -> Ok (`File path)
  in
  Arg.conv ~docv:"File, prefix or directory"
    (fpd_parser check_dir, fpd_printer)

let file_prefix_or_directory_arg action =
  let docv = "FILE|DIRECTORY|PATH-to-PREFIX" in
  let doc =
    let dir_str = sprintf
      "If a directory is specified (eg. path-to-IMGT/alignments/) then by \
      default files from %s loci are added (see --%s). \
      Flags, such as --%s, --%s and --%s determine which files from the loci \
      are selected or merged."
      full_class1_loci_str
      only_class1_flag
      genetic_argument
      nuclear_argument
      merged_argument
    in
    file_prefix_doc ~action ~dir_str
  in
  Arg.(required
      & pos 0 (some fpd_constructor) None
      & info ~doc ~docv [])

let file_prefix_or_directory_to_allele_input_list ?distance ~selectors
  fpd only_class1 gen_nuc_mgd =
  let open Alleles.Input in
  match fpd with
  | `File path   -> [ alignment ?distance ~selectors path]
  | `Prefix p    -> [ merge ?distance ~selectors p]
  | `Directory d ->
      let cls =
        if only_class1 then
          Nomenclature.ClassI
        else
          Nomenclature.FullClassI
      in
      let files = to_files gen_nuc_mgd d cls in
      begin match gen_nuc_mgd with
      | `Merged   -> List.map files ~f:(merge ?distance ~selectors)
      | `Genetic
      | `Nuclear  -> List.map files ~f:(alignment ?distance ~selectors)
      end

let errored error_code fmt =
  ksprintf (fun s -> eprintf "%s" s; error_code) fmt

let at_most_two_fastqs ~single ~paired = function
  | []        -> errored Term.exit_status_cli_error "No FASTQ files specified"
  | [f1]      -> single f1; Term.exit_status_success
  | [f1; f2]  -> paired f1 f2; Term.exit_status_success
  | lst       -> errored Term.exit_status_cli_error
                    "Too many FASTQ files specified: %d"
                    (List.length lst)

let phmm_construction_error = 4
let phmm_construction_exit_info =
  Term.exit_info ~doc:"PHMM construction error" phmm_construction_error

exception Phmm_construction_error of string
