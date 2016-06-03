
let make n e =
  let k = Pattern.pow4 n in
  k, Array.make k e

let init n f =
  let k = Pattern.pow4 n in
  k, Array.init k ~f:(fun i -> f (Pattern.int_to_pat ~k i))

