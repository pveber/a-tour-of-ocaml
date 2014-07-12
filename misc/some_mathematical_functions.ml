(** {2 Simple functions} *)
let square x = x *. x

let abs x = if x > 0. then x else -. x

let arithmetic_mean x y = (x +. y) /. 2.

let geometric_mean x y = sqrt (x *. y)

let eval f x = f x
(* eval has type : ('a -> 'b) -> 'a -> 'b, which means that eval takes
   two arguments, a function from 'a to 'b, a value of type 'a and
   returns a value of type 'b *)

(* alternative equivalent definitions of eval *)
let eval f = fun x -> f x

let eval = fun f -> fun x -> f x

let eval = fun f x -> f x

let compose f g x = f (g x)
(* compose has type ('a -> 'b) -> ('c -> 'a) -> 'c -> 'b, which means
   it takes two functions such that the range of the first is the
   domain of the second, an argument from the domain of the second and
   returns a value from the range of the first. C'est clair, non ? *)

(* alternative strictly equivalent definitions of compose *)

let compose f g = fun x -> f (g x)

let compose f = fun g -> fun x -> f (g x)

let compose = fun f -> fun g -> fun x -> f (g x)

let compose = fun f g x -> f (g x)




(** {2 Functions on 2d vectors} *)

let origin = (0., 0.)
let ex = (1., 0.)
let ey = (0., 1.)

let norm2 (x, y) = sqrt (x ** 2. +. y ** 2.)

let ( +| ) (x1, y1) (x2,y2) = (x1 +. x2, y1 +. y2)

let ( *| ) lambda (x, y) = (lambda *. x, lambda *. y)

let r = norm2 (ex +| 2. *| ey)
(*
  val r : float = 2.23606797749979
*)

(** {2 First steps with recursion} *)
let rec sigma n =
  if n = 0 then 0 else n + sigma (n - 1)
(* this function is unsafe: [sigma (- 1)] never returns (because it
   enters an infinite recursion). In practice, it fails with an
   exception Stack_overflow, which happens when there are too many
   nested function calls. *)

let rec sigma n =
  if n < 0 then failwith "sigma: invalid argument"
  else if n = 0 then 0
  else n + sigma (n - 1)

let rec factorial n =
  if n < 0 then failwith "factorial: invalid argument"
  else if n = 0 then 1
  else n * factorial (n - 1)

let rec power n f =
  if n < 0 then failwith "power: invalid argument"
  else if n = 0 then fun x -> x
  else fun x -> power (n - 1) f (f x)
(* example:

# power 5 (fun x -> x + 1) 0;;
- : int = 5

*)

(** {2 Application: derivatives of a function} *)

let derivative eps f x = (f (x +. eps) -. f x) /. eps
(*
# derivative 1e-6 square 0.;;
- : float = 1e-06
# derivative 1e-6 square 1.;;
- : float = 2.00000099992436731
# derivative 1e-6 square 2.;;
- : float = 4.00000100064801245
*)

let relative_error x ~from:x_ref =
  abs_float (1. -. x /. x_ref)
(* note the use of a labeled argument to distinguish between the
   approximation and the reference value. Both have type float so the
   type check won't help you distinguish them. Using a label, you
   don't need to remember the order of arguments anymore. *)

let numerical_derivative_approximation ~f ~f'' x =
  let numerical_f'' = power 2 (derivative 1e-6) f in
  relative_error (numerical_f'' x) ~from:(f'' x)
(*
# numerical_derivative_approximation ~f:log ~f'':(fun x -> -. 1. /. (x ** 2.)) 1.;;
- : float = 2.0000380231977033e-06
*)

(** {2 Application: computing a root square by dichotomy} *)
let rec dichotomic_square_root_aux eps x a b =
  let m = arithmetic_mean a b in
  if b -. a < eps then m
  else (
    let a', b' = if m *. m > x then a, m else m, b in
    dichotomic_square_root_aux eps x a' b'
  )

let dichotomic_square_root ?(eps = 1e-6) x =
  if eps < 0. || x < 0. then failwith "dichotomic_square_root: invalid argument"
  else dichotomic_square_root_aux eps x 0. x
(*
# dichotomic_square_root 4.;;
- : float = 2.0000004768371582

Note the use of optional arguments. In the previous call the argument
eps is ommited, and is given its default value in the body of
dichotomic_square_root

# dichotomic_square_root ~eps:1e-9 2.;;
- : float = 1.41421356191858649
# sqrt 2.;;
- : float = 1.41421356237309515
# dichotomic_square_root ~eps:1e-9 2. -. sqrt 2.;;;;
- : float = -4.5450865293616971e-10

*)

let rec dichotomic_square_root_aux debug eps x a b =
  let () = if debug then Printf.printf "Current value of the interval: (%f,%f)\n" a b in
  let m = arithmetic_mean a b in
  if b -. a < eps then m
  else (
    let a', b' = if m *. m > x then a, m else m, b in
    dichotomic_square_root_aux debug eps x a' b'
  )

let dichotomic_square_root ?(debug = false) ?(eps = 1e-6) x =
  if eps < 0. || x < 0. then failwith "dichotomic_square_root: invalid argument"
  else dichotomic_square_root_aux debug eps x 0. x

(*
# dichotomic_square_root ~debug:true ~eps:1e-3 3.;;
Current value of the interval: (0.000000,3.000000)
Current value of the interval: (1.500000,3.000000)
Current value of the interval: (1.500000,2.250000)
Current value of the interval: (1.500000,1.875000)
Current value of the interval: (1.687500,1.875000)
Current value of the interval: (1.687500,1.781250)
Current value of the interval: (1.687500,1.734375)
Current value of the interval: (1.710938,1.734375)
Current value of the interval: (1.722656,1.734375)
Current value of the interval: (1.728516,1.734375)
Current value of the interval: (1.731445,1.734375)
Current value of the interval: (1.731445,1.732910)
Current value of the interval: (1.731445,1.732178)
- : float = 1.7318115234375
*)

(** {2 Application: Babylonian method to compute a root square} *)

let rec babylonian_square_root_aux guess n y =
  if n = 0 then guess y
  else
    let x = babylonian_square_root_aux guess (n - 1) y in
    0.5 *. (x +. y /. x)

let babylonian_square_root ?(guess = fun x -> x) n y =
  if n < 0 || y < 0. then failwith "babylonian_square_root: invalid argument"
  else babylonian_square_root_aux guess n y

(*
# sqrt 3.;;
- : float = 1.73205080756887719
# babylonian_square_root 3 3.;;
- : float = 1.73214285714285721
# babylonian_square_root 4 3.;;
- : float = 1.73205081001472738
# babylonian_square_root 5 3.;;
- : float = 1.73205080756887719

The convergence rate of this method is awfully fast: the number of
correct digits doubles at each iteration!

Note that in this implementation, we left the choice of an initial
guess as a parameter (but with a default, trivial value). That way, we
can test several implementations of the initial guess very
easily. Let's do it!
*)
let normalize x =
  let rec loop a n =
    if a < 100. && n mod 2 = 0 then a, n
    else loop (a /. 10.) (n + 1)
  in
  loop x 0
(*
# normalize 125348.;;
- : float * int = (12.5348, 4)
# normalize 12534.;;
- : float * int = (1.2534, 4)
*)

let sqrt_guess x =
  let a, n = normalize x in
  let alpha = if a < 10. then 2. else 6. in
  alpha *. 10. ** (float (n / 2))
(*
# sqrt_guess 125348.;;
- : float = 600.
# sqrt 125348.;;
- : float = 354.045194855120144
# babylonian_square_root 3 125348.;;
- : float = 15671.1249162369295
# babylonian_square_root ~guess:sqrt_guess 3 125348.;;
- : float = 354.059011038189
#

It's indeed much better!
*)

