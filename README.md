# traffic-assignment

Given a network transportation, this code computes two Wardrop Equilibrium, one based on density of cars in each arc of the network (the density is the quantity of 
cars taken in an instant divided by the length of the arc), the other one based on flow of vehicles of cars as classical. The approach based on density allows 
for more realistic cost functions. 

Function [N,NN,D,c,cc,Dcrit,Tw,lambd,mu,err,iter]=asignacioncantidades(q,Delta,Deltaw,alpha,beta,l,L,R,P,DN,Dtilde):

Arguments:
1) q is the matrix which indicates the repartition of the cars in the path of the network. More precisely, if a is an arc which belongs to the path p, then q_{pa} is the number of cars in arc a divided by the number of car in path p. If the arc a does not belong to the path p, then q_{pa}=0. The matrix q can be optimized in a bilevel programming.
2) Delta is the incidence matrix between paths and arcs, that is, Delta_{pa}=1 if the arc a belongs to the path p, Delta_{pa}=0 if not.
3) Deltaw is the incidence matrix between pairs origin-destination and path, that is, Deltaw_{wp}=1 if the path p joins the pair origin-destination w, Deltaw_{wp}=0 if not.
4) alpha and beta are the vector of parameters of the mean speed functions, that is for each arc a, v_a(D_a)=alpha_a - beta_a*D_a, where D_a is the density of cars in arc a. 



Transport1.m , Transport2.m , Transport3.m , Transport4.m , Transport5.m : different configuration of the traffic network.

