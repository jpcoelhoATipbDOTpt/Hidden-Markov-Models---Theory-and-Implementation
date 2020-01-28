function x=rosen(gradF,x,A,b,C,d,maxiter)
% Rosen algorithm for solving:
% min f(X)
% s.t.
%     AX-b<=0
%     CX-d =0

% Initialization phase...
[Aa,Ai,bi,Nk]=MatGradRest(A,b,C,x);

% Gradient calculation of f(x)
df=feval(gradF,x);
for itera=1:maxiter
    disp(['------------------------Iteration: ' num2str(itera)]) 
    if isempty(Nk)
        dk=-df;
        if sum(abs(dk))==0, %dk=0
            % The optimum was found
            break;
        else
            x=update_sol(gradF,x,Ai,b,dk);
            [Aa,Ai,bi,Nk]=MatGradRest(A,b,C,x);
        end
    else
        % Determine the projection matrix...
        %P=Nk'\(Nk*Nk')*Nk;
        P=Nk'*inv(Nk*Nk')*Nk;
        P=P-eye(size(P));
        % Determine search Direction
        dk=P*df;
        % checks whether critical point...
        if sum(abs(dk))==0, % dk=0
            % Calculate the Lagrange multipliers
            lambda=-inv(Nk*Nk')*Nk*df;
            % Checks the signal of the active restrictions associated
            % with inequalities
            if ~isempty(Aa) % If there are active inequality constraints
                [nl,nc]=size(Aa); % nl indicates the number of 
                                  % active constraints
                neglambda=find(lambda(1:nl)<0);
                if ~isempty(neglambda) % there are negative multipliers,
                    [val,ind]=min(lambda);
                    Nk(ind,:)=[]; % Remove the constraint associated with
                                  % the most negative multiplier...
                    Aa(ind,:)=[];
                else % There are no negative multipliers. 
                     % Conditions of KKT are met...
                    break;
                end
            end
        else % dk~=0
            x=update_sol(gradF,x,Ai,bi,dk);
            [Aa,Ai,bi,Nk]=MatGradRest(A,b,C,x);
        end
    end
    df=feval(gradF,x);
end