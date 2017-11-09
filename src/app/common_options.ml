(* Common arguments shared between programs. *)

module Ns = String
open Util
open Cmdliner

let (//) = Filename.concat

let repo = "prohlatype"

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
let alignment_arg =
  let docv = "FILE" in
  let doc  = "File to lookup IMGT allele alignments. The alleles found in this \
              file will initially define the set of alleles to be used. Use an \
              allele selector to modify this set." in
  Arg.(value & opt (some file) None & info ~doc ~docv ["alignment"])

let merge_argument = "merge"

let merge_arg, merges_arg =
  let parser_ path =
    let s = Filename.basename path in
    match Nomenclature.parse_locus s with
    | Error e -> `Error e
    | Ok l ->
        let n = path ^ "_nuc.txt" in
        let g = path ^ "_gen.txt" in
        if not (List.mem ~set:Alter_MSA.supported_loci l) then
          `Error ("gene not supported: " ^ s)
        else if not (Sys.file_exists n) then
          `Error ("Nuclear alignment file doesn't exist: " ^ n)
        else if not (Sys.file_exists g) then
          `Error ("Genetic alignment file doesn't exist: " ^ n)
        else
          `Ok path  (* Return path, and do appending later, the prefix is more useful. *)
  in
  let convrtr = parser_, (fun frmt -> Format.fprintf frmt "%s") in
  let loci_s ~sep =
    List.map ~f:Nomenclature.show_locus Alter_MSA.supported_loci
    |> String.concat ~sep
  in
  let docv = sprintf "[%s]" (loci_s ~sep:"|") in
  let doc  =
    sprintf "Construct a merged (gDNA and cDNA) graph of the specified \
            prefix path. Currently only supports %s genes. The argument must \
            be a path to files with $(docv)_nuc.txt and $(docv)_gen.txt. \
            Combines with the file arguments to determine the set of loci to \
            type at the same time. The set of alleles is defined by the \
            ones in the nuc file."
      (loci_s ~sep:", ")
  in
  Arg.(value & opt (some convrtr) None & info ~doc ~docv ["m"; merge_argument])
  , Arg.(value & opt_all convrtr [] & info ~doc ~docv ["m"; merge_argument])

(*** Allele selector arguments. ***)
let regex_command_line_args = ["allele-regex"]
let allele_command_line_args = ["spec-allele"]
let without_command_line_args = ["without-allele"]
let num_command_line_args = ["n"; "num-alt"]
let do_not_ignore_suffixed_alleles_flags = ["do-not-ignore-suffixed-alleles"]

let args_to_string lst =
  List.map lst ~f:(fun s -> (if String.length s = 1 then "-" else "--") ^ s)
  |> String.concat ~sep:", "

let regex_arg =
  let docv = "REGEX" in
  let doc  = "Specify alleles to add to the graph via a regex. This option is \
              similar to allele, but lets the user specify a wider range \
              (specifically those alleles matching a POSIX regex) of alleles \
              to include. The '*' used in HLA allele names must be properly \
              escaped from the command line: ex. -ar \"A*02:03\". This \
              \"allele selector\" has the highest precedence, and is applied \
              to all of the alleles from the working set"
  in
  let open Arg in
  let parser_ = parser_of_kind_of_string ~kind:docv
    (fun s -> Some (Alleles.Selectors.Regex s))
  in
  value & opt_all (conv ~docv (parser_, Alleles.Selectors.pp)) []
        & info ~doc ~docv regex_command_line_args

let allele_arg =
  let docv = "STRING" in
  let doc  =
    sprintf "Specify specfic alleles to add to the graph. \
             This \"allele selector\" is applied after regex (%s) \
             but before the without (%s). \
             Use this option to construct a graph with specific alleles (ex. \
             A*01:02). One can repeat this argument to specify multiple \
             alleles."
      (args_to_string regex_command_line_args)
      (args_to_string without_command_line_args)
  in
  let open Arg in
  let parser_ =
    parser_of_kind_of_string ~kind:docv
      (fun s -> Some (Alleles.Selectors.Specific s))
  in
  value & opt_all (conv ~docv (parser_, Alleles.Selectors.pp)) []
        & info ~doc ~docv allele_command_line_args

let without_arg =
  let docv = "STRING" in
  let doc  =
    sprintf "Alleles to remove from the working set. \
             This \"allele selector\" is applied before the exclude (%s)
             selector but after specific (%s) one.
             Use this option to construct a graph without alleles (ex. \
             A*01:02). One can repeat this argument to specify multiple \
             alleles to exclude."
      (args_to_string allele_command_line_args)
      (args_to_string num_command_line_args)
  in
  let open Arg in
  let parser_ =
    parser_of_kind_of_string ~kind:docv
      (fun s -> Some (Alleles.Selectors.Without s))
  in
  value & opt_all (conv ~docv (parser_, Alleles.Selectors.pp)) []
        & info ~doc ~docv without_command_line_args

let num_alt_arg =
  let docv = "POSITIVE INTEGER" in
  let doc  =
    sprintf "Number of alternate alleles to add to the graph. \
             If not specified, all of the alternate alleles in the alignment \
             file are added. This \"allele selector\" is applied last to the
             working set of alleles derived from the alignment file \
             (ex. A_gen, DRB_nuc) or merge request (ex. B) after the without
             argument (%s)"
             (args_to_string without_command_line_args)

  in
  let open Arg in
  let parser_ = (positive_int_parser (fun d -> Alleles.Selectors.Number d)) in
  let nconv = conv ~docv (parser_, Alleles.Selectors.pp) in
  (value & opt (some nconv) None & info ~doc ~docv num_command_line_args)

let do_not_ignore_suffixed_alleles_flag =
  let doc = "IMGT publishes alleles with suffixes ('N', 'L', 'S', 'C', 'A' or \
            'Q') to indicate special conditions such as null alleles ('N') or \
            questionable expression ('Q'). Often times, these alleles can \
            complicate prohlatype as reads sometimes align suprisingly well \
            to them (ex. HLA-A*11:50Q) as opposed to the true \
            allele. For the purposes of having a clearer typing algorithm we \
            exclude all such alleles from IMGT by default. For rare cases or \
            diagnostic purposes one can pass this flag to bring them back \
            into the analysis."
  in
  Arg.(value & flag & info ~doc do_not_ignore_suffixed_alleles_flags)

(*** Other args. ***)
(*
let remove_reference_flag =
  let doc  = "Remove the reference allele from the graph. The reference \
              allele is the one that is listed first in the alignments file. \
              Graphs are currently constructed based upon their \"diff\" to \
              the reference as represented in the alignments file. Therefore, \
              the original reference sequence must be a part of the graph \
              during construction. Specifying this flag will remove it from \
              the graph after the other alleles are added."
  in
  Arg.(value & flag & info ~doc ["no-reference"])
  *)

let no_cache_flag =
  let doc =
    sprintf "Do not use a disk cache (in %s sub directory of the current \
             directory) to search for previously (and then save) constructed \
             graphs."
      Cache.dir
  in
  Arg.(value & flag & info ~doc ["no-cache"])

let do_not_join_same_sequence_paths_flag =
  let doc = "Do not join same sequence paths; remove overlapping nodes at the \
             same position, in the string graph."
  in
  Arg.(value & flag & info ~doc ["do-not-join-same-sequence-paths"])

let option_to_list o =
  Option.value_map o ~default:[] ~f:(fun s -> [s])

let aggregate_selectors ?number_alleles
  ~regex_list ~specific_list ~without_list ~do_not_ignore_suffixed_alleles =
    regex_list
    @ specific_list
    @ without_list
    @ (option_to_list number_alleles)
    @ (if do_not_ignore_suffixed_alleles then
        [Alleles.Selectors.DoNotIgnoreSuffixed] else [])

let to_allele_input ?alignment_file ?merge_file ?distance ~selectors =
  let open Alleles.Input in
  match merge_file with
  | Some prefix   -> Ok (merge ?distance ~selectors prefix)
  | None  ->
      match alignment_file with
      | Some path -> Ok (alignment ?distance ~selectors path)
      | None      -> Error "Neither a file nor merge argument was specified."

(* The more than one case. *)
let to_allele_inputs ~alignment_files ~merge_files ?distance ~selectors =
  let als =
    List.map alignment_files ~f:(Alleles.Input.alignment ?distance ~selectors)
  in
  let mls =
    List.map merge_files ~f:(Alleles.Input.merge ?distance ~selectors)
  in
  match als @ mls with
  | []     -> Error "Neither a merge nor alignment file specified!"
  | inputs -> Ok inputs

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
  Arg.(non_empty & pos_all file [] & info ~doc ~docv [])

let num_reads_arg =
  let docv = "POSITIVE INTEGER" in
  let doc = "Number of reads to take from the front of the FASTA file" in
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
  let doc  = "Read name string (to be found in fastq) to type. Add multiple \
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
              seen this many mismatches." in
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

let map_depth_argument = "map-depth"
let map_depth_default = 5

let allele_depth_arg =
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Specify a positive integer to indicate the number of best alleles \
             and positions to report. Defaults to %d. Be cautious about \
             increasing this number as keeping track of this data can slow \
             down the final analysis."
      map_depth_default
  in
  Arg.(value & opt positive_int map_depth_default
             & info ~doc ~docv [map_depth_argument])

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
  let doc  =
    sprintf "Tweak the accuracy (and consequently speed) of the forward pass \
      (larger -> less accurate but faster). The default value is %f. The \
      implementation uses 63 bit floats which is more than sufficient for \
      this calculation."
      !ParPHMM.dx
  in
  Arg.(value & opt (some positive_float) None
             & info ~doc ~docv ["numerical-accuracy"])

let do_not_past_threshold_filter_flag =
  let doc = "Do not use previously calculated likelihoods (ex. from other loci \
             or considering the reverse complement alignment) as a threshold \
             filter to short-circuit evaluations."
  in
  Arg.(value & flag & info ~doc ["do-not-past-threshold-filter"])

let likelihood_first_flag =
  let doc = "When printing the final per-allele likelihoods, print ordered \
              by likelihood first. This will print all of the alleles, in a \
              compressed format following the likelihood."
  in
  Arg.(value & flag & info ~doc ["likelihood-first"; "llhdfst"])

let do_not_finish_singles_flag =
  let doc = "Paired FASTQ files that may be created by downstream analysis may \
             not guarantee that all of the reads are paired. When in paired \
             mode (2 FASTQ files are passed), by default we'll continue \
             evaluating the unpaired reads as if they were eingle. By passing \
             this flag we'll skip them entirely and use only the paired reads."
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
  let open ParPHMM_drivers.Zygosity_array in
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
        | NonZero v -> Format.fprintf frmt "%f" v)

let zygosity_report_size_argument = "zygosity-report-size"
let zygosity_report_size_arg =
  let open ParPHMM_drivers in
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Override the default number of allelic pairs reported as part of \
      the zygosity portion. By defaults only the pairs that have non-zero \
      (> %f) probability will be reported. Set this value to zero to not \
      report any values."
      Zygosity_array.default_non_zero
  in
  Arg.(value & opt non_negative_zygosity
                      Zygosity_array.(NonZero default_non_zero)
             & info ~doc ~docv [zygosity_report_size_argument])

let zygosity_non_zero_value_arg =
  let open ParPHMM_drivers in
  let docv = "POSITIVE FLOAT" in
  let doc =
    sprintf "Override the default lowerbound of non-zero zygosities. The \
      default (%0.5f) probability might be too high for some scenarios, such \
      as if there are too few reads. This argument overrides %s."
      Zygosity_array.default_non_zero zygosity_report_size_argument
  in
  Arg.(value & opt (some positive_float) None
             & info ~doc ~docv ["zygosity-non-zero"])

let to_num_zygosities ~zygosity_non_zero_value ~zygosity_report_size =
  Option.value_map zygosity_non_zero_value
    ~f:(fun v -> ParPHMM_drivers.Zygosity_array.NonZero v)
    ~default:zygosity_report_size

let per_reads_report_size_arg =
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Override the number of (per)-read information to report. By \
      default all of the per-read information is reported."
  in
  Arg.(value & opt (some non_negative_int) None
             & info ~doc ~docv ["per-read-report-size"])

let number_processes_arg =
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Parallelize read processing across this number of processesors. \
             By default the application is not parallelized as the user would \
             most likely achieve greater efficiency by natively parallelizing \
             across samples as opposed to parallelizing across reads. \
             Furthermore the user MUST specify %S as well in this mode."
      read_size_override_argument
  in
  Arg.(value & opt (some positive_int) None & info ~doc ~docv ["number-processors"])

let split_arg =
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Split the processing of a read into this many smaller forward \
      passes, that are then aggregated. This should decrease the run time cost \
      of a full pass as the match probabilities, across alleles, are less \
      distributed. This value $(b, must) divide the read\ size evenly. \
      Becareful about splitting too much (setting this value too high) as \
      longer reads align more accurately than shorter reads (specifically, \
      10bp is too short while 25bp rarely gets filtered). If the resulting \
      output results in too many filtered reads, consider using a smaller \
      value."
  in
  Arg.(value & opt (some positive_int) None & info ~doc ~docv ["split"])

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
  let doc = "Change output from stdout to files that are prefixed with the \
    specified argument. A log file will be generated with the the \".log\" \
    suffix. Data output will be written to a file with either the \".tsv\" \
    or \".json\" suffix."
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

let class1gen_arg =
  let docv  = "DIRECTORY" in
  let doc   = "Short-cut argument that expands the given dir to look for \
                A_gen.txt, B_gen.txt and C_gen.txt. This overrides any \
                alignment files or merge arguments." in
  Arg.(value & opt (some dir) None & info ~doc ~docv ["class1-gen"])

let class1nuc_arg =
  let docv  = "DIRECTORY" in
  let doc   = "Short-cut argument that expands the given dir to look for \
                A_gen.txt, B_gen.txt and C_gen.txt. This overrides any \
                alignment files or merge arguments." in
  Arg.(value & opt (some dir) None & info ~doc ~docv ["class1-nuc"])

let class1mgd_arg =
  let docv  = "DIRECTORY" in
  let doc   = "Short-cut argument that expands the given dir to look for \
                A_gen.txt, B_gen.txt and C_gen.txt. This overrides any \
                alignment files or merge arguments." in
  Arg.(value & opt (some dir) None & info ~doc ~docv ["class1-mgd"])

(* These are convenience selectors that choose more than one locus at a time. *)
let class_selectors class1_gen_dir class1_nuc_dir class1_mgd_dir
  alignment_files merge_files =
  match class1_gen_dir, class1_nuc_dir, class1_mgd_dir with
  | Some d, _,      _       ->
      (d // "A_gen.txt") :: (d // "B_gen.txt") :: (d // "C_gen.txt") :: alignment_files
      , merge_files
  | None,   Some d, _       ->
      (d // "A_nuc.txt") :: (d // "B_nuc.txt") :: (d // "C_nuc.txt") :: alignment_files
      , merge_files
  | None,   None,   Some d  ->
      alignment_files
      , (d // "A") :: (d // "B") :: (d // "C") :: merge_files
  | None,   None,   None    ->
      alignment_files
      , merge_files
