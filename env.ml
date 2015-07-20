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
   Average ages of death assuming no other possible cause of death
*)
let mean_age_of_death_a = 100.
let mean_age_of_death_a_plus = 200.

let mean_age_of_death_b = 100.
let mean_age_of_death_b_plus = 200.

let mean_age_of_death_other = 40.

let lambda_a = 1. /. mean_age_of_death_a
let lambda_a_plus = 1. /. mean_age_of_death_a_plus

let lambda_b = 1. /. mean_age_of_death_b
let lambda_b_plus = 1. /. mean_age_of_death_b_plus

let lambda_other = 1. /. mean_age_of_death_other

type genome = {
  a_plus: bool;
  b_plus: bool;
}

let string_of_genome x =
  sprintf "{%s %s}"
    (if x.a_plus then "A+" else "A ")
    (if x.b_plus then "B+" else "B ")

let print_report gen cause age =
  printf "%s cause %s, age %.2f\n"
    (string_of_genome gen) cause age

let sim_one gen =
  (*
     Obtain age of an individual when first lethal accident
     for each cause of death (A, B, other) occurs
  *)
  let a =
    exponential_distribution
      (if gen.a_plus then lambda_a_plus else lambda_a)
  in
  let b =
    exponential_distribution
      (if gen.b_plus then lambda_b_plus else lambda_b)
  in
  let other =
    exponential_distribution lambda_other
  in
  let cause_of_death, age =
    let sorted =
      List.sort (fun (_, x) (_, y) -> compare x y)
        ["A", a;
         "B", b;
         "other", other]
    in
    match sorted with
    | (cause, age) :: _ -> cause, age
    | _ -> assert false
  in
  print_report gen cause_of_death age
