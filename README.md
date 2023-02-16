# traffic-assignment

Given a network transportation, this code computes two Wardrop Equilibrium, one based on density of cars in each arc of the network (the density is the quantity of 
cars taken in an instant divided by the length of the arc), the other one based on flow of vehicles of cars as classical. The approach based on density allows 
for more realistic cost functions. 

Function [N,NN,D,c,cc,Dcrit,Tw,lambd,mu,err,iter]=asignacioncantidades(q,Delta,Deltaw,alpha,beta,l,L,R,P,DN,Dtilde):

Arguments of the function:
1) q is the matrix which indicates the repartition of the cars in the path of the network. More precisely, if a is an arc which belongs to the path p, then q_{pa} is the number of cars in arc a divided by the number of car in path p. If the arc a does not belong to the path p, then q_{pa}=0. The matrix q can be optimized in a bilevel programming.
2) Delta is the incidence matrix between paths and arcs, that is, Delta_{pa}=1 if the arc a belongs to the path p, Delta_{pa}=0 if not.
3) Deltaw is the incidence matrix between pairs origin-destination and path, that is, Deltaw_{wp}=1 if the path p joins the pair origin-destination w, Deltaw_{wp}=0 if not.
4) alpha and beta are the vector of parameters of the mean speed functions, that is for each arc a, v_a(D_a)=alpha_a - beta_a*D_a, where D_a is the density of cars in arc a. 
5) l is the vector of the lengths of the arcs.
6) L, R and P are respectively the number of arcs, pairs origin-destination and paths in the network.
7) DN is the vectors of numbers of cars travelling between each pair origin-destination.

Image of the function:
1) N is the vector of numbers of cars travelling in each path.
2) NN is the vector of numbers of cars travelling in each arc. We have the relation N = q*NN.
3) D is the density of cars in each arc, D_a=NN_a/l_a.
4) c is the vector of the time travel in each path.
5) cc is the vector of the time travel in each arc. We have c=Delta*cc.
6) Dcrit is the density which maximizes the flow in each arc, it is given by Dcrit_a = alpha_a/(2*beta_a).
7) lambd is the Lagrange multiplier associated with the satisfaction of the demand, it can be interpreted as a marginal cost with respect to the number of cars travelling between each pair origin-destination.
8) mu is the Lagrange multiplier associated with the nonnegativity of the number of cars travelling on each path. It can be interpreted as an extra cost for travelling the unused path (mu_p=0 if N_p>0).
9) err is the norm of the gap function associated to the Wardrop Equilibrium.
10) iter is the number of iterations.

Function [F,FF,cF,ccF,lamb,mu,err]=asignacionflujos(Delta,Deltaw,alpha,beta,l,L,R,P,Tw,Hcost):

The arguments are the same than previous function except the following:
1) Tw is the vector of flows between the pairs origin-destinations.
3) Hcost 



Transport1.m , Transport2.m , Transport3.m , Transport4.m , Transport5.m : different configuration of the traffic network.

