disp("PART 1: Data Normalized by Variance")
disp(DataVar)
V = table2array(DataVar) %normalized data

disp("Transpose of V")
Vtrans = transpose(V) %take transpose

disp("Number of Players")
N = 17 %number of players

disp("Covariance Matrix")
Sv=(1/(N-1))*Vtrans*V %7x7 covariance matrix for data normalized by variance

disp("Eigenvalues and Eigenvectorss of Covariance Matrix")
evalSv = eig(Sv); %eigenvectors for m
[P,D]=eig(Sv); %eigenvector matrix for Sv

%organize diagnolized matrix and eigenvectors in descending order
[d,ind] = sort(diag(D), "descend");
Ds = D(ind,ind)
Ps = P(:,ind)

totVariance = sum(evalSv(:,1)); %total variance(sum of eigenvectors of S)
W= ['Total Variance =',totVariance];
disp(W)
%check total variance
T=trace(Sv);
disp('The trace of matrix S is equal to the total variance which is equal to the sum of the eigenvalues which is equal to:')
disp(T)%find component variances

%making eigenvectors unit eigenvectors
mag = sqrt(sum(Ps'.^2))';
A=[];
for i= 1:numel(mag)
    A(i,:) = Ps(i,:)./mag(i);
end

  disp("Unit eigenvectors:")
disp(A)
 
disp("Principal component(unit eigenvector), followed by corresponding eigenvalue and it's component variance(%)")


for i=1:7
CompVariance = round(Ds(i,i)/totVariance*100,3);
eval =Ds(i,i);
cvariance = CompVariance;
principalcomponent = A(:,i)
X =sprintf('eigenvalue = %s, component variance = %d %' , eval, cvariance) ;
disp(X)
end



%% 


disp("PART 2: Data Normalized by Unit Magnitude")
disp(DataMag)


M = table2array(DataMag)

disp("Transpose of M")
Mtrans = transpose(M)

disp("Number of Players")
N = 17 %number of players

disp("Covariance Matrix")
Sm = (1/(N-1))*Mtrans*M %7x7 covariance matrix for data normalized by magnitude

disp("Eigenvalues and Eigenvectorss of Covariance Matrix")
evalSm = eig(Sm); %eigenvectors for m
[P,D]=eig(Sm); %eigenvector matrix for Sm

%organize diagnolized matrix and eigenvectors in descending order
[d,ind] = sort(diag(D), "descend");
Ds = D(ind,ind)
Ps = P(:,ind)


totVariance = sum(evalSm(:,1)) %total variance(sum of eigenvectors of S)
L= ["Total Variance =",totVariance];
disp(L)
%check total variance
T=trace(Sm);
disp('The trace of matrix S is equal to the total variance which is equal to the sum of the eigenvalues which is equal to:')
disp(T)%find component variances

%making eigenvectors unit eigenvectors
mag = sqrt(sum(Ps'.^2))';
A=[];
for i= 1:numel(mag)
    A(i,:) = Ps(i,:)./mag(i);
end

disp("Unit eigenvectors:")
disp(A)

%find component variances
disp("Principal component(unit eigenvector), followed by corresponding eigenvalue and it's component variance(%)")

for i=1:7
CompVariance = round(Ds(i,i)/totVariance*100,3);
eval =Ds(i,i);
cvariance = CompVariance;
principalcomponent = A(:,i)
X =sprintf('eigenvalue = %s,  component variance = %d' , eval, cvariance) ;
disp(X)
end

