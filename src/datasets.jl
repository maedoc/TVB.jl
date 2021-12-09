using DelimitedFiles
using ZipFile


function conn76()
    fname = "conn76.zip"
    if ~isfile(fname)
        url = "https://github.com/the-virtual-brain/tvb-data/raw/master/tvb_data/connectivity/connectivity_76.zip"
        download(url, fname)
    end
    zf = ZipFile.Reader(fname)
    # @show zf.files
    l = readdlm(zf.files[6])
    w = readdlm(zf.files[7])
    close(zf)
    w, l
end