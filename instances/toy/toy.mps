NAME ./toy
OBJSENSE 
 MAX
ROWS
 N cost
 G f0
 G f1
 G f2
 G f3
COLUMNS
 x cost -1
 x f0 25
 x f1 -1
 x f2 -2
 x f3 2
 y cost -10
 y f0 -20
 y f1 -2
 y f2 1
 y f3 10
RHS
 rhs f0 -30
 rhs f1 -10
 rhs f2 -15
 rhs f3 15
BOUNDS
 LI  bound  x	4.00
 UI  bound  x	8.00
 LB  bound  y	0.00
ENDATA

