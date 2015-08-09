(*
   Test whether it is faster to select an allele of one gene
   and then an allele of another gene via successive environments
   or whether a single environment that favors both is faster or
   just as fast.

   gene A, alleles a+ and a-
   gene B, alleles b+ and b-

   environment A: favors a+, neutral for other genes
   environment B: favors b+, neutral for other genes
   environment AB: favors both a+ and b+, neutral for other genes

   Assumptions:

   - genes considered influence the lethality of accidents that can
     happen at any time during the individual's life. The probability
     of such lethal accidents is modeled by an exponential distribution.

   - reproduction is asexual; parent eventually dies of old age.
     (uh, we need to find a way to recombine genomes)

   - the number of children is proportional to the age of death. The idea
     is that general health results both in more children and longer life.

   - even though individuals die at different ages, we assume that the
     age difference between parent and children (one generation) is constant.

   - the population is kept constant by normalizing the number of children
     from generation to generation.

   Examples:

   - gene A: influences running speed, critical to survive bear attacks
   - gene B: influences resistance against flu infections
*)

open Printf

let random_bool f =
  Random.float 1. < f

let uniform_distribution () =
  Random.float 1.

(*
   The exponential distribution is continuous and monotonous, so can
   simulate it by taking the inverse of the cumulative distribution function:
   X = F^(-1)(U)   where F(x) = 1 - e^(-lambda * x)
*)
let exponential_distribution lambda =
  -. log (uniform_distribution ()) /. lambda

(*
   Average ages of death for one gene
   assuming no other possible cause of death
*)
let mean_age_of_death_minus = 100.
let mean_age_of_death_plus = 110.

let lambda_minus = 1. /. mean_age_of_death_minus
let lambda_plus = 1. /. mean_age_of_death_plus

type genome = bool array

let sob = function
  | true -> "1"
  | false -> "0"

let string_of_genome x =
  sprintf "{%s}"
    (String.concat "" (List.map sob (Array.to_list x)))

let print_individual_report gen i age =
  printf "%s cause %i(%s), age %.2f\n"
    (string_of_genome gen)
    i (sob gen.(i))
    age

let print_population_report gens =
  assert (Array.length gens > 0);
  let gene_count = Array.length gens.(0) in
  let plus_fractions =
    Array.init gene_count (fun i ->
      let n =
        Array.fold_left (fun acc gen ->
          if gen.(i) then acc + 1
          else acc
        ) 0 gens
      in
      float n /. float (Array.length gens)
    )
  in
  Array.iteri (fun i x ->
    printf "gene %i: %.2f%%\n" i (100. *. x)
  ) plus_fractions

let mini a =
  assert (Array.length a > 0);
  BatArray.fold_lefti (fun ((_, m) as best) i v ->
    if v < m then (i, v)
    else best
  ) (0, a.(0)) a

let sim_one gen =
  (*
     Obtain age of an individual when first lethal accident
     for each cause of death (one of the genes) occurs
  *)
  let ages =
    Array.mapi (fun i g ->
      let lambda =
        match g with
        | true -> lambda_plus
        | false -> lambda_minus
      in
      let age = exponential_distribution lambda in
      age
    ) gen
  in
  let i, age_of_death = mini ages in
  print_individual_report gen i age_of_death;
  age_of_death

(* Scale elements of the array such that their sum equals 1 *)
let normalize a =
  let sum = Array.fold_left (+.) 0. a in
  Array.map (fun x -> x /. sum) a

(* In-place swapping of some fraction of the genes between
   random individuals *)
let recombine gens0 =
  let gens = Array.map Array.copy gens0 in
  let pop_size = Array.length gens in
  assert (pop_size > 0);
  let gene_count = Array.length gens.(0) in
  for i = 0 to pop_size - 1 do
    for g = 0 to gene_count - 1 do
      if random_bool 0.1 then
        let j = Random.int pop_size in
        if i <> j then
          let tmp = gens.(i).(g) in
          gens.(i).(g) <- gens.(j).(g);
          gens.(j).(g) <- tmp
    done
  done;
  gens

(* Simulate a population over max_generations *)
let rec sim_pop max_generations generation target_pop_size gens =
  printf "--- Generation %i ---\n%!" generation;
  printf "Population size: %i\n" (Array.length gens);
  print_population_report gens;
  let ages = Array.map sim_one gens in
  let relative_weights = normalize ages in
  let child_counts = Array.map (fun w ->
    truncate (w *. float target_pop_size +. 0.5)
  ) relative_weights
  in
  let new_pop_size = Array.fold_left (+) 0 child_counts in
  assert (new_pop_size >= 1);
  let cloned_children =
    Array.concat (
      Array.to_list (
        Array.mapi (fun i gen ->
          Array.make child_counts.(i) (Array.copy gen)
        ) gens
      )
    )
  in
  let children = recombine cloned_children in
  if generation < max_generations then
    sim_pop max_generations (generation + 1) target_pop_size children

let main () =
  let max_generations = 1000 in
  let gene_count = 10 in
  let plus_allele_freq = 0.1 in
  let target_pop_size = 1000 in
  let gens =
    Array.init target_pop_size (fun _ ->
      Array.init gene_count (fun _ -> random_bool plus_allele_freq)
    )
  in
  sim_pop max_generations 1 target_pop_size gens

let () = main ()
