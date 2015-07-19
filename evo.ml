(*
  Multi-population evolution, with constant fitness function within
  each population.
*)

open Printf

type config = {
  num_genes : int; (* number of genes *)
  high_weight : int; (* weight of major gene in the fitness function *)
  low_weight : int; (* weight of other genes in the fitness function *)

  pop_size : int; (* number of individuals in a population *)
  recomb_ratio : float; (* probability of a gene being swapped with another
                           individual of the same population *)

  num_pop : int; (* number of populations, each having a different major 
                    gene *)
  migration_ratio : float; (* probability of an individual to be swapped
                              with an individual from another population *)
}

type gene_instance = int (* 0 or 1 *)

type individual = gene_instance array

type population = {
  fitness_weights : int array;
  individuals : individual array;
}


let init_pop conf i =
  let num_indiv = conf.pop_size in
  let num_genes = conf.num_genes in
  let major_gene = i in
  if major_gene >= conf.num_genes then
    failwith "Too many populations or not enough genes for the model";
  let genes = 
    Array.init num_genes
      (fun i ->
         if i = major_gene then conf.high_weight
         else conf.low_weight
      )
  in
  let a =
    Array.init num_indiv (
      fun i ->
        Array.init num_genes (
          fun j ->
            (* initial frequency of the better allele is 10% *)
            if Random.int 10 = 0 then 1
            else 0
        )
    )
  in
  {
    fitness_weights = genes;
    individuals = a;
  }

let init_species conf =
  Array.init conf.num_pop (init_pop conf)


let deaths conf p = 
  let a = p.individuals in
  let num_indiv = Array.length a in
  let b = Array.init conf.pop_size (fun _ -> a.(Random.int num_indiv)) in
  { p with individuals = b }

let fitness weights x =
  let f = ref 0 in
  for i = 0 to Array.length weights - 1 do
    f := !f + weights.(i) * x.(i)
  done;
  !f

let recombine conf p =
  let a = p.individuals in
  let num_indiv = Array.length a in
  let num_genes = Array.length p.fitness_weights in
  let num_recomb =
    truncate (0.5 *. conf.recomb_ratio *. float (num_indiv * num_genes))
  in
  (* swap 2 random genes from 2 random individuals *)
  for k = 1 to num_recomb do
    let x = a.(Random.int num_indiv) in
    let y = a.(Random.int num_indiv) in
    let g = Random.int num_genes in
    let xg = x.(g) in
    x.(g) <- y.(g);
    y.(g) <- xg
  done;
  p
  

let births p =
  let l =
    Array.fold_left (
      fun accu x ->
        let f = fitness p.fitness_weights x in
        (Array.make f x) :: accu
    ) [] p.individuals
  in
  { p with individuals = Array.concat l }


let migrate conf sp =
  let num_pop = conf.num_pop in
  let num_indiv = conf.pop_size in
  let num_swap =
    truncate (0.5 *. conf.migration_ratio *. float (num_pop * num_indiv))
  in
  (* swap 2 random individuals from 2 random populations *)
  for k = 1 to num_swap do
    let x = sp.(Random.int num_pop).individuals in
    let y = sp.(Random.int num_pop).individuals in
    let i = Random.int num_indiv in
    let j = Random.int num_indiv in
    let xi = x.(i) in
    x.(i) <- y.(j);
    y.(j) <- xi
  done;
  sp


let cycle conf sp =
  for i = 0 to Array.length sp - 1 do
    sp.(i) <- recombine conf (deaths conf (births sp.(i)));
    assert (Array.length sp.(i).individuals = conf.pop_size);
  done;
  migrate conf sp


let frequencies p =
  let num_genes = Array.length p.fitness_weights in
  let a = Array.make num_genes 0 in
  Array.iter (
    fun x ->
      for i = 0 to Array.length x - 1 do
        a.(i) <- a.(i) + x.(i)
      done
  ) p.individuals;
  let num_indiv = Array.length p.individuals in
  a, num_indiv

let average_frequency counts n =
  Array.fold_left (+) 0 counts, n * Array.length counts


let print_frequencies a n =
  let nf = float n in
  for i = 0 to Array.length a - 1 do
    printf "[%i] %g\n" i (float a.(i) /. nf)
  done

let add_array accu a =
  assert (Array.length accu = Array.length a);
  for i = 0 to Array.length accu - 1 do
    accu.(i) <- accu.(i) + a.(i)
  done



let stop_condition counts nmax stop =
  try 
    for i = 0 to stop - 1 do
      if counts.(i) < nmax then
        raise Exit
    done;
    true
  with Exit -> false


let print_cycle_info conf i sp =
  printf "--- Cycle %i ---\n" i;
  let num_pop = Array.length sp in
  let num_genes = Array.length sp.(0).fitness_weights in
  let global_counts, global_n = Array.make num_genes 0, ref 0 in
  for j = 0 to num_pop - 1 do
    let p = sp.(j) in
    printf "population %i: %i\n" j (Array.length p.individuals);
    let counts, counts_n = frequencies p in
    add_array global_counts counts;
    global_n := !global_n + counts_n;
    let av_count, av_count_n = average_frequency counts counts_n in
    printf "frequency of major gene %i: %g\n"
      j (float counts.(j) /. float counts_n);
    printf "average frequency: %g\n"
      (float av_count /. float av_count_n);
    (*print_frequencies freq;*)
  done;
  (*printf "global frequencies:\n";
  print_frequencies global_counts !global_n;*)
  flush stdout;
  stop_condition global_counts !global_n conf.num_pop



let run conf =
  let sp0 = init_species conf in
  let sp = ref sp0 in
  let cycle_id = ref 0 in
  try
    let stop = print_cycle_info conf !cycle_id !sp in
    if stop then
      raise Exit;
    while true do
      incr cycle_id;
      sp := cycle conf !sp;
      let stop = print_cycle_info conf !cycle_id !sp in
      if stop then
        raise Exit
    done
  with Exit ->
    printf "Converged and succeeded after %i generations.\n" !cycle_id


let conf = {
  num_genes = 10;
  high_weight = 100;
  low_weight = 1;
  pop_size = 1000;
  recomb_ratio = 0.5;
  num_pop = 10;
  migration_ratio = 0.01;
}

let () =
  run conf
