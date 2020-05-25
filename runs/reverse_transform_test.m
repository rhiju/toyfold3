ctr1 = ctr(:,1);
ctr2 = ctr(:,2);
M1 = M(:,:,1);
M2 = M(:,:,2);

[tf,Rf] = get_transform( ctr1, M1, ctr2, M2);

[tr,Rr] = get_transform( ctr2, M2, ctr1, M1)

[tr_rev, Rr_rev] =  reverse_transform( tf,Rf)