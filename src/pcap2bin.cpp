// Convert a (gzipped) pcap trace into the binary stream format, one record per IPv4 packet.
//   pcap2bin <trace.pcap[.gz]> <out.bin> [max_packets]
// Attributes: srcIP, dstIP, srcPort, dstPort, proto, ipLen, src/16, dst/16, ttl, tcpFlags
#include "common.hpp"

static inline uint16_t be16(const unsigned char* b) { return (uint16_t)((b[0] << 8) | b[1]); }
static inline uint32_t be32(const unsigned char* b) { return ((uint32_t)b[0] << 24) | ((uint32_t)b[1] << 16) | ((uint32_t)b[2] << 8) | b[3]; }

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: pcap2bin trace out [max]\n");
        return 1;
    }
    std::string path = argv[1];
    uint64_t maxp = argc > 3 ? strtoull(argv[3], 0, 10) : UINT64_MAX;
    bool gz = path.size() > 3 && path.substr(path.size() - 3) == ".gz";
    FILE* f = gz ? popen(("gzip -dc '" + path + "'").c_str(), "r") : fopen(path.c_str(), "rb");
    if (!f) return 1;
    unsigned char gh[24];
    if (fread(gh, 1, 24, f) != 24) return 1;
    uint32_t magic;
    memcpy(&magic, gh, 4);
    bool swap = (magic == 0xd4c3b2a1);
    if (magic != 0xa1b2c3d4 && !swap) {
        fprintf(stderr, "not a classic pcap file (magic %08x)\n", magic);
        return 1;
    }
    auto rd32 = [&](const unsigned char* b) {
        uint32_t x;
        memcpy(&x, b, 4);
        return swap ? __builtin_bswap32(x) : x;
    };
    uint32_t linktype = rd32(gh + 20);
    fprintf(stderr, "linktype=%u\n", linktype);
    Dataset d;
    d.p = 10;
    std::vector<unsigned char> buf(1 << 16);
    uint64_t total = 0, kept = 0;
    unsigned char ph[16];
    while (kept < maxp && fread(ph, 1, 16, f) == 16) {
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
        } else if (linktype != 101 && linktype != 12) {
            continue;
        }
        if (incl < off + 20 || (b[off] >> 4) != 4) continue;
        const unsigned char* ip = b + off;
        uint32_t ihl = (ip[0] & 15u) * 4;
        uint32_t proto = ip[9], src = be32(ip + 12), dst = be32(ip + 16);
        uint32_t sport = 0, dport = 0, flags = 0;
        if ((proto == 6 || proto == 17) && incl >= off + ihl + 4) {
            sport = be16(ip + ihl);
            dport = be16(ip + ihl + 2);
            if (proto == 6 && incl >= off + ihl + 14) flags = ip[ihl + 13];
        }
        uint32_t rec[10] = {src, dst, sport, dport, proto, be16(ip + 2), src >> 16, dst >> 16, ip[8], flags};
        d.v.insert(d.v.end(), rec, rec + 10);
        kept++;
    }
    d.n = kept;
    d.save(argv[2]);
    fprintf(stderr, "packets=%llu ipv4_kept=%llu\n", (unsigned long long)total, (unsigned long long)kept);
    return 0;
}
