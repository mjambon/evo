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
