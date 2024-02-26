function [contourX,contourY]=sort_by_nearest(contourX,contourY);
    %for complicated shapes, follow the points one by one
    points_left=1;
    newX=[contourX(1)];
    newY=[contourY(1)];
    X_todo=contourX(2:end);
    Y_todo=contourY(2:end);
    while points_left
        curX=newX(end);
        curY=newY(end);
        axs=1:length(X_todo);
        [~, idx]=min(((X_todo-curX).^2+(Y_todo-curY).^2).^0.5);
         idx=idx(1);
         newX(length(newX)+1)=X_todo(idx);
         newY(length(newY)+1)=Y_todo(idx);
         X_todo=X_todo(axs~=idx);
         Y_todo=Y_todo(axs~=idx);
         if isempty(X_todo)
             points_left=0;
         end
    end
    contourX=newX';
    contourY=newY';