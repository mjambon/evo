(*
  Single-population evolution, with fitness function varying
  over time.
*)

open Printf

type config = {
  num_genes : int; (* number of genes *)
  high_weight : float; (* weight of major gene in the fitness function *)
  low_weight : float; (* weight of other genes in the fitness function *)

  pop_size : int; (* number of individuals in the population *)
  recomb_ratio : float; (* probability of a gene being swapped with another
                           individual of the same population *)
  births_per_cycle : float; (* average number of births per cycle
                               per individual *)

  period : int; (* period during which the major gene remains the same *)
}

type gene_instance = Good | Bad

type individual = gene_instance array

type population = {
  major_gene : int;
  fitness_weights : float array;
  individuals : individual array;
}

let make_fitness_weights conf g =
  Array.init conf.num_genes 
    (fun i ->
       if i = g then conf.high_weight
       else conf.low_weight)


let update_fitness_weights conf pop i =
  if i mod conf.period = 0 then
    let g = (i / conf.period) mod conf.num_genes in
    { pop with
        major_gene = g;
        fitness_weights = make_fitness_weights conf g }
  else
    pop

let init_pop conf =
  let num_indiv = conf.pop_size in
  let num_genes = conf.num_genes in
  let weights = make_fitness_weights conf 0 in
  let a =
    Array.init num_indiv (
      fun i ->
        Array.init num_genes (
          fun j ->
            (* initial frequency of the better allele is 10% *)
            if Random.int 10 = 0 then Good
            else Bad
        )
    )
  in
  {
    major_gene = 0;
    fitness_weights = weights;
    individuals = a;
  }

(* pick n random elements from array a (note: elements of a are permuted) *)
let pick a n =
  if n > Array.length a then
    invalid_arg "pick";
  let m = Array.length a - n in
  let b = Array.make n a.(0) in
  for i = n - 1 downto 0 do
    let j = Random.int (m+i+1) in
    b.(i) <- a.(j);
    a.(j) <- a.(m+i);
  done;
  b

let deaths conf p = 
  let a = p.individuals in
  let num_indiv = Array.length a in
  if num_indiv <= conf.pop_size then p
  else
    { p with individuals = pick a conf.pop_size }

let fitness weights x =
  let f = ref 0. in
  for i = 0 to Array.length weights - 1 do
    if x.(i) = Good then
      f := !f +. weights.(i)
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
  

let births conf pop =
  let raw_fitness_values =
    Array.map (fun x -> fitness pop.fitness_weights x) pop.individuals in
  let raw_fitness_sum =
    Array.fold_left (+.) 0. raw_fitness_values in
  let normalizer =
    conf.births_per_cycle *. float (Array.length pop.individuals)
    /. raw_fitness_sum
  in
  let normalized_fitness_values =
    Array.map (fun x -> normalizer *. x) raw_fitness_values
  in
  let children =
    let accu = ref [] in
    for i = 0 to Array.length normalized_fitness_values - 1 do
      let x = normalized_fitness_values.(i) in
      (* x = expected number of children *)
      let nmin = truncate x in
      let p = x -. float nmin in
      let extra_child = Random.float 1. < p in
      let n =
        if extra_child then nmin + 1
        else nmin in
      if n > 0 then
        accu := Array.make n pop.individuals.(i) :: !accu
    done;
    Array.concat !accu
  in
  let expected_births =
    conf.births_per_cycle *. float (Array.length pop.individuals) in
  let actual_births = Array.length children in
  printf "%i births (expected %g)\n" actual_births expected_births;
  { pop with individuals = (*Array.append pop.individuals*) children }


let cycle conf p i =
  let p = recombine conf (deaths conf (births conf p)) in
  update_fitness_weights conf p i


let frequencies p =
  let num_genes = Array.length p.fitness_weights in
  let a = Array.make num_genes 0 in
  Array.iter (
    fun x ->
      for i = 0 to Array.length x - 1 do
        match x.(i) with
            Good -> a.(i) <- a.(i) + 1
          | Bad -> ()
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



let stop_condition counts nmax =
  try 
    for i = 0 to Array.length counts - 1 do
      if counts.(i) > 0 && counts.(i) < nmax then
        raise Exit
    done;
    true
  with Exit -> false


let print_cycle_info conf p i =
  printf "--- Cycle %i ---\n" i;
  printf "population: %i\n" (Array.length p.individuals);
  let counts, counts_n = frequencies p in
  let av_count, av_count_n = average_frequency counts counts_n in
  let major = p.major_gene in
  printf "frequency of major gene %i: %g\n"
    major (float counts.(major) /. float counts_n);
  printf "average frequency: %g\n"
    (float av_count /. float av_count_n);
  (*print_frequencies counts counts_n;*)
  flush stdout;
  stop_condition counts counts_n



let run conf =
  let p0 = init_pop conf in
  let p = ref p0 in
  let cycle_id = ref 0 in
  try
    let stop = print_cycle_info conf !p !cycle_id in
    if stop then
      raise Exit;
    while true do
      incr cycle_id;
      p := cycle conf !p !cycle_id;
      let stop = print_cycle_info conf !p !cycle_id in
      if stop then
        raise Exit
    done
  with Exit ->
    printf "Converged after %i generations.\n" !cycle_id


let conf = {
  num_genes = 100;
  high_weight = 100.;
  low_weight = 1.;
  pop_size = 1000;
  recomb_ratio = 0.1;
  births_per_cycle = 1.1;
  period = 20;
}

let () =
  Random.self_init ();
  run conf
