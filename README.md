# calc
Command-line mathematical expression and equation evaluator

## How to use
To run either supply an expression or equation as a command line parameter
```
./calc "4/3 * pi * 7.5^3"
= 1767.15
```
Or run without arguments to continuously enter expressions
```
./calc
> 4/3 * pi * 7.5^3
= 1767.15
> x^2 - 5x + 6 = 0
x = 3
> exit
```

## Supported functions/constants
Function/constant | token
--- | ---
Add | `+`
Subtract | `-`
Multiply | `*`
Power | `^`
Sine | `sin()`
Cosine | `cos()`
Tangent | `tan()`
Inverse Sine | `asin()`
Inverse Cosine | `scos()`
Inverse Tangent | `atan()`
Square Root | `sqrt()`
Absolute | `abs()`
Logarithm | `log()`
Natural Logarithm | `ln()`
Pi | `pi`
e | `e`
Previous Answer | `ans`
