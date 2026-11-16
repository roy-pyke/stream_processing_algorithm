// Convert (gzipped) pcap traces into the binary stream format, one record per IPv4 packet.
//   pcap2bin <trace.pcap[.gz]> <out.bin> [max_packets]            MAWI: 10 attributes
//       srcIP, dstIP, srcPort, dstPort, proto, ipLen, src/16, dst/16, ttl, tcpFlags
//   pcap2bin --caida <out.bin> <trace.pcap[.gz]>...                CAIDA: 6 attributes, traces in the given order
//       dstIP, proto, dst/16, srcIP, src&0xffff, src/16
// The CAIDA attributes are those of the OmniSketch paper (omniVLDB.CompareBaselines, attribute set
// {ipDst, ipProto, ipDstNet, ipSrc, ipSrcHost, ipSrcNet}); the authors encode an address by deleting its
// dots, which merges some addresses (1.23.4.5 and 12.3.4.5), whereas we keep the exact 32-bit address
// and the exact halves. Link types: Ethernet (1, VLAN tags skipped) and raw IPv4 (12, 101, 228);
// microsecond and nanosecond pcap files of either byte order.
#include "common.hpp"

static inline uint16_t be16(const unsigned char* b) { return (uint16_t)((b[0] << 8) | b[1]); }
static inline uint32_t be32(const unsigned char* b) { return ((uint32_t)b[0] << 24) | ((uint32_t)b[1] << 16) | ((uint32_t)b[2] << 8) | b[3]; }

// Appends the IPv4 packets of one trace to d; returns false if the file cannot be read as pcap.
static bool convert(const std::string& path, bool caida, Dataset& d, uint64_t maxp, uint64_t& total) {
    bool gz = path.size() > 3 && path.substr(path.size() - 3) == ".gz";
    FILE* f = gz ? popen(("gzip -dc '" + path + "'").c_str(), "r") : fopen(path.c_str(), "rb");
    if (!f) return false;
    unsigned char gh[24];
    if (fread(gh, 1, 24, f) != 24) return false;
    uint32_t magic;
    memcpy(&magic, gh, 4);
    bool swap = (magic == 0xd4c3b2a1 || magic == 0x4d3cb2a1);
    if (magic != 0xa1b2c3d4 && magic != 0xa1b23c4d && !swap) {
        fprintf(stderr, "%s: not a classic pcap file (magic %08x)\n", path.c_str(), magic);
        return false;
    }
    auto rd32 = [&](const unsigned char* b) {
        uint32_t x;
        memcpy(&x, b, 4);
        return swap ? __builtin_bswap32(x) : x;
    };
    uint32_t linktype = rd32(gh + 20);
    std::vector<unsigned char> buf(1 << 16);
    unsigned char ph[16];
    uint64_t kept0 = d.n;
    while (d.n < maxp && fread(ph, 1, 16, f) == 16) {
        uint32_t incl = rd32(ph + 8);
        if (incl > buf.size()) buf.resize(incl);
        if (fread(buf.data(), 1, incl, f) != incl) break;
        total++;
        const unsigned char* b = buf.data();
        uint32_t off = 0;
        if (linktype == 1) { // Ethernet, possibly VLAN-tagged
            if (incl < 14) continue;
            uint16_t et = be16(b + 12);
            off = 14;
            while (et == 0x8100 && incl >= off + 4) et = be16(b + off + 2), off += 4;
            if (et != 0x0800) continue;
        } else if (linktype != 101 && linktype != 12 && linktype != 228) {
            continue;
        }
        if (incl < off + 20 || (b[off] >> 4) != 4) continue;
        const unsigned char* ip = b + off;
        uint32_t ihl = (ip[0] & 15u) * 4;
        uint32_t proto = ip[9], src = be32(ip + 12), dst = be32(ip + 16);
        if (caida) {
            uint32_t rec[6] = {dst, proto, dst >> 16, src, src & 0xffffu, src >> 16};
            d.v.insert(d.v.end(), rec, rec + 6);
        } else {
            uint32_t sport = 0, dport = 0, flags = 0;
            if ((proto == 6 || proto == 17) && incl >= off + ihl + 4) {
                sport = be16(ip + ihl);
                dport = be16(ip + ihl + 2);
                if (proto == 6 && incl >= off + ihl + 14) flags = ip[ihl + 13];
            }
            uint32_t rec[10] = {src, dst, sport, dport, proto, be16(ip + 2), src >> 16, dst >> 16, ip[8], flags};
            d.v.insert(d.v.end(), rec, rec + 10);
        }
        d.n++;
    }
    if (gz) pclose(f); else fclose(f);
    fprintf(stderr, "%s: linktype=%u ipv4_kept=%llu\n", path.c_str(), linktype, (unsigned long long)(d.n - kept0));
    return true;
}

int main(int argc, char** argv) {
    bool caida = argc > 1 && std::string(argv[1]) == "--caida";
    if (argc < 3 || (caida && argc < 4)) {
        fprintf(stderr, "usage: pcap2bin trace out [max]  |  pcap2bin --caida out trace...\n");
        return 1;
    }
    Dataset d;
    d.p = caida ? 6 : 10;
    uint64_t total = 0;
    std::string out;
    if (caida) {
        out = argv[2];
        for (int i = 3; i < argc; i++)
            if (!convert(argv[i], true, d, UINT64_MAX, total)) return 1;
    } else {
        out = argv[2];
        if (!convert(argv[1], false, d, argc > 3 ? strtoull(argv[3], 0, 10) : UINT64_MAX, total)) return 1;
    }
    d.save(out);
    fprintf(stderr, "packets=%llu ipv4_kept=%llu\n", (unsigned long long)total, (unsigned long long)d.n);
    return 0;
}
