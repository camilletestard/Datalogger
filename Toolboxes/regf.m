function errors = regf(X1train,ytrain,X1test,ytest)
mdl = fitlm(X1train, ytrain);
yfit = predict(mdl,X1test);
MAE = mean(abs(yfit-ytest));
adjMAE = MAE/range(ytest);
Rsq = corr(yfit,ytest)^2;
errors = [MAE adjMAE Rsq];
end