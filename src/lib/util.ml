
include Nonstd
module String = Sosa.Native_string

let invalid_argf fmt = ksprintf invalid_arg fmt

let disk_memoize ?dir arg_to_string f =
  let dir = Option.value dir ~default:(Filename.get_temp_dir_name ()) in
  fun ?(skip_disk_cache=false) arg ->
    if skip_disk_cache then
      f arg
    else
      let file = Filename.concat dir (arg_to_string arg) in
      if Sys.file_exists file then begin
        let i = open_in file in
        let r = Marshal.from_channel i in
        close_in i;
        r
      end else begin
        if Sys.command (sprintf "mkdir -p %s" dir) <> 0 then
          invalid_argf "Failed to make sure cache dir %s is available" dir;
        let r = f arg in
        let o = open_out file in
        Marshal.to_channel o r [];
        close_out o;
        r
      end

let error_bind ~f oe =
  match oe with
  | Error e -> Error e
  | Ok o -> f o

let (>>=) oe f = error_bind ~f oe

let error_map ~f oe =
  match oe with
  | Error e -> Error e
  | Ok o    -> Ok (f o)

let (>>|) oe ~f = error_map ~f oe

let error fmt = ksprintf (fun s -> Error s) fmt

let unwrap_ok = function
  | Ok o    -> o
  | Error _ -> invalid_arg "Not Ok in unwrap_ok."

let unwrap_error = function
  | Ok _    -> invalid_argf "Not Error in unwrap_error."
  | Error e -> e

let short_seq s =
  let n = String.length s in
  if n > 10 then
    sprintf "%s...%s"
      (String.slice_exn ~finish:4 s) (String.slice_exn ~start:(n-3) s)
  else
    s

let index_string s index =
  let (b, a) = String.split_at s ~index in
  b ^ "." ^ a

(** Compare two strings and display bars for vertical mismatches. *)
let manual_comp_display s1 s2 =
  let msm = ref 0 in
  let cs = String.mapi s1 ~f:(fun index c1 ->
    let c2 = String.get_exn s2 ~index in
    if c1 = c2 then ' ' else
      begin incr msm; '|' end)
  in
  let ms = string_of_int !msm in
  let n  = String.length ms + 1 in
  let pd = String.make n ' ' in
  sprintf "%s%s\n%s %s\n%s%s"
    pd s1 ms cs pd s2
