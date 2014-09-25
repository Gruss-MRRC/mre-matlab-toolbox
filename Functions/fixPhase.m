function im_out = fixPhase(im_in,p1,sl,fr)

thresh= 300;

s1= size(im_in);


im_out = im_in;
for j = 1:s1(1),
    for k = 1:p1,
        im_intens = im_in(j,k,sl,fr,1);
        del_im = im_in(j,k+1,sl,fr,2)-im_in(j,k,sl,fr,2);
        if (im_intens > thresh && del_im < -2048),
            im_out(j,k+1:end,sl,fr,2) = im_out(j,k+1:end,sl,fr,2) + 4096;
        end
    end
    for k = p1+1:s1(2)-1,
        im_intens = im_in(j,k,sl,fr,1);
        del_im = im_in(j,k+1,sl,fr,2)-im_in(j,k,sl,fr,2);
        if (im_intens > thresh && del_im > 2048),
            im_out(j,k+1:end,sl,fr,2) = im_out(j,k+1:end,sl,fr,2) - 4096;
        end
    end
end
            