
val align_sequences :
      s1:string ->
      o1:int ->
      s2:string ->
      o2:int -> [> `Both of int | `First of int * int | `Second of int * int ]

type 'e pa = { mismatches : int; edge : 'e; alignment : 'e alignment; }
and 'e alignment = Leaf of int | Partial of 'e pa list

val map_alignment : ('a -> 'b) -> 'a alignment -> 'b alignment
val align :
      ?mub:int ->
      Ref_graph.G.t ->
      Graph_index.t ->
      string -> (Ref_graph.Edges.t alignment list, string) result

val name_edges_in_alignment :
      Alleles.index -> Alleles.Set.t alignment -> string alignment
val best_of_paths : 'a alignment -> int
val to_weights : int list -> float list
val init_alingment_map : Alleles.index -> float Alleles.Map.t
val alignments_to_weights :
      ?base_weight:float ->
      float Alleles.Map.t -> Alleles.Set.t alignment -> unit
val most_likely :
      Alleles.index -> float Alleles.Map.t -> (float * string) list
