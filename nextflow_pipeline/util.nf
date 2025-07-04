
def header(label, offset=8) {
    def off = (1..offset).collect{ ' ' }.join()
    def hdr = label.replaceAll(/./,'=')
    return "${off}${hdr}\n${off}${label}\n${off}${hdr}"
}