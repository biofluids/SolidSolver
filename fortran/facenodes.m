function list=facenodes(nsd,nen,face)
list = zeros(no_facenodes(nsd,nen),1);
i3 = [2,3,1];
i4 = [2,3,4,1];
if nsd==2
    if nen==3
        list = [face;i3(face)];
    elseif nen==4
        list = [face;i4(face)];
    end
elseif nsd==3
    if nen==4
        if face==1
            list = [1,2,3];
        elseif face==2
            list = [1,4,2];
        elseif face==3
            list = [2,4,3];
        elseif face==4
            list = [3,4,1];
        end
    elseif nen==8
        if face==1
            list = [1,2,3,4];
        elseif face==2
            list = [5,8,7,6];
        elseif face==3
            list = [1,5,6,2];
        elseif face==4
            list = [2,6,7,3];
        elseif face==5
            list = [3,7,8,4];
        elseif face==6
            list = [4,8,5,1];
        end
    end
end
end