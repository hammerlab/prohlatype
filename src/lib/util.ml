
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

