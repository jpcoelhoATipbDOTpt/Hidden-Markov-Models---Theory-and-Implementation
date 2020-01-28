function [Aa,Ai,bi,Nk,x]=one_shot_rosen(gradfunc,Theta,O,A,b,C,d,Aa,Ai,bi,Nk,x,df,var)

if isempty(Nk)
    dk=-df;
    if sum(abs(dk))==0, %dk=0
        % The optimum was found
        return;
    else
        x=hmmupdate_sol(gradfunc,Theta,O,x,Ai,bi,dk,var);
        [Aa,Ai,bi,Nk]=MatGradRest(A,b,C,x);
    end
else
    % Determine the projection matrix...
    P=Nk'*((Nk*Nk')\Nk);
    P=P-eye(size(P));
    % Determine search direction
    dk=P*df;
    % checks whether critical point...
    if sum(abs(dk))==0, % dk=0
        % Calculate the Lagrange multipliers
        lambda=-inv(Nk*Nk')*Nk*df;
        % Checks the signal of the active restrictions associated with inequalities
        if ~isempty(Aa) % If there are active inequality constraints
            [nl,~]=size(Aa); % nl indicates the number of active restrictions
            neglambda=find(lambda(1:nl)<0);
            if ~isempty(neglambda) % there are negative multipliers...
                [~,ind]=min(lambda);
                Nk(ind,:)=[]; % remove the constraint associated with the most negative multiplier
                Aa(ind,:)=[];
            else % there are no negative multipliers. KKT conditions are met...
                return;
            end
        end
    else % dk~=0
        x=hmmupdate_sol(gradfunc,Theta,O,x,Ai,bi,dk,var);
        [Aa,Ai,bi,Nk]=MatGradRest(A,b,C,x);
    end
end