(** Computing with probabilities is a pain, we'll encode different Rings to
    hide the actual implementation. *)

open Util

(*let dx = 2.22044e-16 *)
let dx = ref 1.e-6

let is_nan x = x <> x

let close_enough x y =
  if is_nan x then
    if is_nan y then
      true
    else
      false
  else if is_nan y then
    false
  else
    let d = x -. y in
    (abs_float d) < !dx

(* Probability Ring where we perform the forward pass calculation. *)
module type Ring = sig

  type t

  val pp : Format.formatter -> t -> unit
  val to_yojson : t -> Yojson.Safe.json
  val of_yojson : Yojson.Safe.json -> (t, string) result
  val to_string : ?precision:int -> t -> string
  val zero : t
  val one  : t

  val gap  : t
  val is_gap : t -> bool

  val ( + ) : t -> t -> t
  val ( * ) : t -> t -> t
  val ( / ) : t -> t -> t
  val max   : t -> t -> t
  val ( < ) : t -> t -> bool
  val ( <= ) : t -> t -> bool
  val compare : t -> t -> int

  val close_enough : t -> t -> bool
  val close_enough_cmp : t -> t -> int

  (* Special constructs necessary for the probabilistic logic. *)
  (* Convert constant probabilities. *)
  val constant : float -> t

  (* Scale a probability be a third. *)
  val times_one_third : float -> t

  (* Complement probability. *)
  val complement_probability : float -> t

  val probability : ?maxl:t -> t -> float

  val as_float : ?precision:int -> t -> float

  val pow : float -> t -> t

end (* Ring *)

module Just_floats = struct

  (* Code that is common to both rings. *)
  let gap   = nan
  let is_gap = is_nan

  let ( < ) x y =
    if is_gap x then true else
      if is_gap y then false else
        x < y

  let ( <= ) x y =
    if is_gap x then true else
      if is_gap y then false else
        x <= y

  let max x y =
    if x <= y then y else x

  let compare x y =
    if is_gap x then -1 else
      if is_gap y then 1 else
        compare x y

  (* For some weird reason this isn't inlined correctly and it is much faster
   * to leave the original close_enough defined above. *)
  let close_enough x y =
    close_enough x y

  let close_enough_cmp x y =
    if close_enough x y then
      0
    else if x < y then
      -1
    else
      1

  let to_string ?(precision=10) t =
    sprintf "%.*f" precision t

  let as_float ?precision x =
    Option.value_map precision ~default:x
      ~f:(fun p ->
            let m = 10. ** (float p) in
            (floor (x *. m)) /. m)

end (* Just_floats *)

(* Mostly for demonstration purposes. *)
module RegularFloats = struct

  type t = float [@@deriving show,yojson]

  include Just_floats

  let zero  = 0.
  let one   = 1.

  let ( + ) = ( +. )
  let ( * ) = ( *. )
  let ( / ) = ( /. )

  let constant x = x

  let complement_probability p =
    1. -. p

  let times_one_third p =
    p /. 3.

  let probability ?maxl x = x

  let pow x t = t ** x

end (* RegularFloats *)

(* For the calculations inside of ParPHMM we'll just use Log10 because the
   errors are already represented that way via Phred quality scores. A
   benchmark left TODO is to benchmark this choice. *)
module Log10 : Ring = struct

  type t = float
  [@@deriving show,yojson]

  include Just_floats

  let zero  = neg_infinity
  let one   = 0.  (* log10 1. *)

  let gap   = nan
  let is_gap = is_nan

  let exp10 x = 10. ** x

  let ( * ) lx ly = lx +. ly
  let ( / ) lx ly = lx -. ly

  let ( + ) lx ly =
         if lx = neg_infinity then ly
    else if ly = neg_infinity then lx
    else if lx > ly           then lx +. log10 (1. +. exp10 (ly -. lx))
    else (* lx < ly *)             ly +. log10 (1. +. exp10 (lx -. ly))

  let max   = max

  let constant = log10

  let l13 = constant (1. /. 3.)

  let times_one_third = ( * ) l13

  let complement_probability lq =
    log10 (1. -. (exp10 lq))

  (* Embed the log-sum-exp procedure in converting back to probabilities *)
  let probability ?maxl l =
    match maxl with
    | None    -> exp10 l
    | Some ml -> exp10 (l -. ml)

  let pow x t = x *. t

  (* The base error (qualities) are generally know. To avoid repeating the manual
    calculation (as described above) of the log quality to log (1. -. base error)
    we precompute these values.

  Weird. this seems to be slower! TODO: Why? )
  let log10_one_minus_l = function
    | -0.0 -> log_z
    | -0.1 -> -0.686825324380115454
    | -0.2 -> -0.432923433336248276
    | -0.3 -> -0.302062439928300397
    | -0.4 -> -0.220480830541908562
    | -0.5 -> -0.165088538626769726
    | -0.6 -> -0.125627577491815079
    | -0.7 -> -0.0966528953262047186
    | -0.8 -> -0.07494036743261491
    | -0.9 -> -0.058435173882679825
    | -1.0 -> -0.0457574905606751153
    | -1.1 -> -0.035944514242268806
    | -1.2 -> -0.0283047837831196247
    | -1.3 -> -0.0223306727357915694
    | -1.4 -> -0.0176431456736382448
    | -1.5 -> -0.0139554338820558448
    | -1.6 -> -0.0110483332892353306
    | -1.7 -> -0.00875292940206854816
    | -1.8 -> -0.0069382318574496421
    | -1.9 -> -0.00550215071190342936
    | -2.0 -> -0.00436480540245008826
    | -2.1 -> -0.00346349774554599588
    | -2.2 -> -0.00274889425384098815
    | -2.3 -> -0.00218210128532180967
    | -2.4 -> -0.0017324081870171721
    | -2.5 -> -0.00137553579921727916
    | -2.6 -> -0.00109227082153636758
    | -2.7 -> -0.000867397043708781281
    | -2.8 -> -0.000688856394105166097
    | -2.9 -> -0.000547088803770739681
    | -3.0 -> -0.000434511774017691684
    | -3.1 -> -0.000345109452404739883
    | -3.2 -> -0.000274107777278245887
    | -3.3 -> -0.000217717413117151276
    | -3.4 -> -0.000172930172032690271
    | -3.5 -> -0.000137357693108739246
    | -3.6 -> -0.000109103544996691523
    | -3.7 -> -8.66617872714958135e-05
    | -3.8 -> -6.88364918576594357e-05
    | -3.9 -> -5.46778777877239756e-05
    | -4.0 -> -4.34316198075056039e-05
    | -4.1 -> -3.44986070951026853e-05
    | -4.2 -> -2.7402993817532012e-0
    |    x -> (*eprintf "asked to compute log10_one_minus_l of %f\n" x; *)
              log10_one_minus_l_manual x
 *)

end (* Log10 *)

(* Note that functorizing this code, with

module type Exponentiable = sig
  val exp : float -> float
  val log : float -> float
end

and

module MakeLog(E : Exponentiable) : Ring = struct
   ...

end (* MakeLog *)

module Ln = MakeLog(struct
   let exp = exp
   let log = log
   end)

Seems to incur a big ~8% overhead that the compiler doesn't get rid of
automatically (haven't tested with flambda). It is unfortunate, but in this
case we will Repeat Ourselves. This is the reason that we expliclty name
exp10 in the previous module and DON'T in this case.
*)

module Ln : Ring = struct

  type t = float
  [@@deriving show,yojson]

  include Just_floats

  let zero  = neg_infinity
  let one   = 0.  (* log 1. *)

  let gap   = nan
  let is_gap = is_nan

  let ( * ) lx ly = lx +. ly
  let ( / ) lx ly = lx -. ly

  let ( + ) lx ly =
         if lx = neg_infinity then ly
    else if ly = neg_infinity then lx
    else if lx > ly           then lx +. log (1. +. exp (ly -. lx))
    else (* lx < ly *)             ly +. log (1. +. exp (lx -. ly))

  let max   = max

  let constant = log

  let l13 = constant (1. /. 3.)

  let times_one_third = ( * ) l13

  let complement_probability lq =
    log (1. -. (exp lq))

  (* Embed the log-sum-exp procedure in converting back to probabilities *)
  let probability ?maxl l =
    match maxl with
    | None    -> exp l
    | Some ml -> exp (l -. ml)

  let pow x t = x *. t

end (* Ln *)
