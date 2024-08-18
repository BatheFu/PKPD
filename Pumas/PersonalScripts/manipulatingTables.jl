using DataFrames, DataFramesMeta, PharmaDatasets

df = dataset("demographics_1")
first(df, 5)

df[1, 3]
df[5, "SCR"]

df[10:20, ["ID", "WEIGHT"]]

using Statistics
describe(df)

describe(df, :median, :nunique, cols=["ID", "AGE"])

function my_stat(col)
    count(>(50), col)
end

describe(
    df, :median, :mean, my_stat => "gt50"
)

my_df = @orderby(df, :eGFR)
first(my_df, 5)
## the parentheses can be omitted 
my_df = @orderby df :eGFR
first(my_df, 5)

# Or the function form

@orderby df begin
    :ISMALE
    :eGFR
end

first(my_df, 5)

"""
Use a - before the column name to sort in reverse
"""
my_df = @orderby df begin
    :ISMALE
    -:eGFR
end

first(my_df, 5)

#select

@select df :ID :AGE

#need to wrap them in $() to make them work

propertynames(df)

(@select df $(Not(:ID))) |> propertynames

(@select df $(Between(:SCR, :eGFR))) |> propertynames

(@select df $(r"E")) |> propertynames #with regular expressions

(@select df begin
    $(Cols(:ID, r"E"))
    $(Between(:SCR, :eGFR))
end) |> propertynames

#Broadcasting
my_df = @select df begin
    :ID
    :AGE_in_month = :AGE .* 12
end

using CairoMakie
hist(df.eGFR)

function zscore(v)
    μ = mean(v)
    σ = std(v)
    return [(x - μ) / σ for x in v]
end

df_z = @select df :ID :eGFR_z = zscore(:eGFR)
first(df_z, 5)

hist(df_z.eGFR_z)

# allow us to substitute complex expressions using begin blocks.
my_df = @select df begin
    :ID 
    #whitespace usually doesnt matter 
    :eGFR_z = begin
        μ = mean(:eGFR)
        σ = std(:eGFR)
        return [(x - μ) / σ for x in :eGFR]
    end
end

my_df = @transform df begin
    :eGFR_z = begin
        μ = mean(:eGFR)
        σ = std(:eGFR)
        [(x - μ) / σ for x in :eGFR]
    end
end
first(my_df, 5)
"""
transform is better as it keeps all original columns and values
before you modify.
"""

my_df = @transform df begin
    :eGFR_z = begin
        μ = mean(:eGFR)
        σ = std(:eGFR)
        [(x - μ) / σ for x in :eGFR]
    end
    @byrow :AGE_in_months = :AGE * 12
end
first(my_df, 5)

my_df = @rename df :serum_creatinine = :SCR
#:new = :old 
first(my_df, 5)

#with whitespace
my_df = @rename df $("Serum Creatinine") = :SCR
first(my_df, 5)

#accept functions also
my_df = rename(lowercase, df)
first(my_df, 2)

function customRename1(x)
    replace(x, "A" => "X")
end

my_df = rename(customRename1, df)
first(my_df, 2) # A to X 

#subset; it is row operation
my_df = @subset df begin
    :SCR .> mean(:SCR)
end

first(my_df, 5)

my_df = @rsubset(df, :SCR > 1.2)
first(my_df, 5)

my_df = @subset df begin
    :SCR .> mean(:SCR)
    :eGFR .< median(:eGFR)
end;
first(my_df, 5)

#advanced operations
my_df = @rtransform df @astable begin
    AGE_log = log(:AGE)
    WEIGHT_lbs = :WEIGHT * 2.2
    :AGE_log_WEIGHT_lbs = AGE_log + WEIGHT_lbs
end
first(my_df, 5)

#operations within groups
gdf = groupby(df, :ISMALE)
my_df = wz =  @transform gdf begin
    :WEIGHT_z = zscore(:WEIGHT)
end

first(my_df, 5)

fig = Figure()
not_tr = Axis(fig[1, 1], title = "Not transformed", xlabel = "Weight")
tr = Axis(fig[1, 2], title = "Z-score", xlabel = "Z-score")
w = hist!(not_tr, @rsubset(wz, :ISMALE == 0).WEIGHT)
m = hist!(not_tr, @rsubset(wz, :ISMALE == 1).WEIGHT)
hist!(tr, @rsubset(wz, :ISMALE == 0).WEIGHT_z)
hist!(tr, @rsubset(wz, :ISMALE == 1).WEIGHT_z)
Legend(
    fig[2, 1:2],
    [w, m],
    ["women", "men"],
    orientation = :horizontal,
    tellwidth = false,
    tellheight = true,
)
fig

# a chain version
my_df = @chain df begin
    groupby(:ISMALE)
    @transform :WEIGHT_z = zscore(:WEIGHT)
end
first(my_df, 5)

#combine
my_df = @combine gdf begin
    :AGE_μ = mean(:AGE)
    :WEIGHT_μ = mean(:WEIGHT)
    :total = length(:ID)
    :high_eGFR = count(>(80), :eGFR)
end
first(my_df, 5)

#chain version combine
my_df = @chain df begin
    groupby(:ISMALE)
    @combine begin
        :AGE_μ = mean(:AGE)
        :WEIGHT_μ = mean(:WEIGHT)
        :total = length(:ID)
        :high_eGFR = count(>(80), :eGFR)
    end
end
first(my_df, 5)

#Above is equal to 
my_df = @by df :ISMALE :AGE_μ = mean(:AGE) :WEIGHT_μ = mean(:WEIGHT)
first(my_df, 5)

# IS equal to : 
my_df = @by df  :ISMALE begin
    :AGE_μ = mean(:AGE)
    :WEIGHT_μ = mean(:WEIGHT)
    :total = length(:ID)
    :high_eGFR = count(>(80), :eGFR)
end
first(my_df,2 )