# Functions for numerical statistics
using Dates
using Statistics
using DataFrames

function get_confusion(h,v, hmask, vmask)
    h_binary = h .> .5
    v_binary = v .> .5
    trues = [hmask,vmask]
    masks = [h_binary, v_binary]
    numerical_stats = zeros((2,5))
    for k in 1:length(trues)
        est = masks[k]
        true_mask = trues[k]
        mask_height, mask_width = size(est)
        zero_as_one = 0
        one_as_zero = 0
        one_as_one = 0
        zero_as_zero = 0
        total = 0
        for i in 1:mask_height
            for j in 1:mask_width
                if true_mask[i,j] == 0 && est[i,j] == 0
                    zero_as_zero += 1
                elseif true_mask[i,j] == 0 && est[i,j] == 1
                    zero_as_one += 1
                elseif true_mask[i,j] == 1 && est[i,j] == 0
                    one_as_zero += 1
                else
                    one_as_one += 1
                end
                total += 1
            end
        end
        #print(zero_as_one, " ", one_as_zero, " ", one_as_one, " ", zero_as_zero, ", ", total, "\n")
        numerical_stats[k, 1:5] = [zero_as_zero, one_as_one, zero_as_one, one_as_zero, total]

    end
    return numerical_stats
end


function get_f_score(confusion)
    "gets overall f score for h and v put together. This is the one number that
    describing the binary classifiers performance.
    Recall: proportion of those classified as edge that are actually edges.
    Precision: proportion of all edges classified as edges"
    d1 = sum(confusion[:, 2] + confusion[:, 3])
    if d1 > 0
        recall = sum(confusion[:, 2])/sum(confusion[:, 2] + confusion[:, 3])
    else
        recall = 0
    end
    precision = sum(confusion[:, 2])/sum(confusion[:, 2] + confusion[:, 4])
    if precision + recall > 0
        return 2*(precision * recall)/(precision + recall)
    else
        return 0
    end
end




# ### Test out Functions
# include("objective3.jl")
# img = load("./images2/im5.png");
# imgg = Gray.(img);
# im = convert(Array{Float64}, imgg);
#
# img_hmask = load("./images2/im5_hmask.png");
# imgg_hmask = Gray.(img_hmask);
# im_hmask = convert(Array{Float64}, imgg_hmask);
#
# img_vmask = load("./images2/im5_vmask.png");
# imgg_vmask = Gray.(img_vmask);
# im_vmask = convert(Array{Float64}, imgg_vmask);
#
# imheight, imwidth = size(im)
#
# # ###### Starting Point of optimization
# f = .5 .+ randn(imheight, imwidth)/10;
# h = .5 .+ randn(imheight-1, imwidth)/10;
# v = .5 .+ randn(imheight, imwidth-1)/10;
#
# # ### Interior Point
# niters = 7
# ρs = [1.5625*.2^k for k in 2:niters]
# αws = repeat([.0005], niters-1) # .005 seems to work with p, .002 seems to work well with q. .05 maximum
# # alphas try .0005 to .01
# λs = repeat([.5], niters - 1)
# f,h,v = interior_point(gg, ∇gg, q, ∇q, f, h, v, λs, ρs, αws, im)
#
# confusion = get_confusion(h,v, im_hmask, im_vmask)
# #print(confusion, "\n")
# print("f_score = ", get_f_score(confusion))
#
# mosaic(colorview(Gray, im),
# colorview(Gray, f),
# colorview(Gray, v .> .5),
# colorview(Gray, im_vmask),
# colorview(Gray, h .> .5),
# colorview(Gray, im_hmask))
