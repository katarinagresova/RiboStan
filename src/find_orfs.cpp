#include <Rcpp.h>
using namespace Rcpp;

// Find all ORFs in a single nucleotide string.
//
// Scans `seq` for every occurrence of `start_codon`, then walks in-frame
// (3 nt at a time) until a stop codon is encountered.  Results are appended
// into the caller-supplied vectors `starts`, `ends`, and `indices`.
//
// Positions are 1-based and include both the start and stop codons.
// `min_body` is the minimum number of nucleotides between the end of the start
// codon and the beginning of the stop codon (0 = start + stop adjacent).
static void scan_seq(const std::string& seq,
                     const std::string& start_codon,
                     const std::vector<std::string>& stop_codons,
                     int min_body,
                     int seq_index,
                     std::vector<int>& starts,
                     std::vector<int>& ends,
                     std::vector<int>& indices) {
  int n = (int)seq.size();
  int sc_len = (int)start_codon.size();   // always 3
  std::size_t pos = 0;

  while ((pos = seq.find(start_codon, pos)) != std::string::npos) {
    int s = (int)pos;                       // 0-based start of start codon
    int body_start = s + sc_len;            // 0-based start of codon after ATG

    // Walk in-frame looking for a stop codon
    int p = body_start;
    while (p + 3 <= n) {
      const std::string codon = seq.substr(p, 3);
      bool is_stop = false;
      for (const auto& sc : stop_codons) {
        if (codon == sc) { is_stop = true; break; }
      }
      if (is_stop) {
        int body_len = p - body_start;      // nt between start and stop
        if (body_len >= min_body) {
          starts.push_back(s + 1);          // convert to 1-based
          ends.push_back(p + 3);            // 1-based inclusive end
          indices.push_back(seq_index);
        }
        break;
      }
      p += 3;
    }
    ++pos;  // advance past current position to find all (possibly overlapping) starts
  }
}


//' Find ORFs in a character vector of nucleotide sequences
//'
//' @param seqs Character vector of uppercase DNA sequences.
//' @param start_codon Length-1 character; start codon (default "ATG").
//' @param stop_codons Character vector; stop codons.
//' @param min_body Integer; minimum body length in nucleotides between the end
//'   of the start codon and the beginning of the stop codon.
//' @return A list with three integer vectors: \code{starts}, \code{ends}
//'   (1-based, inclusive), and \code{indices} (1-based index into \code{seqs}).
//' @keywords internal
// [[Rcpp::export]]
List find_orfs_cpp(CharacterVector seqs,
                   std::string     start_codon,
                   CharacterVector stop_codons,
                   int             min_body = 0) {
  // Convert stop_codons to std::vector<std::string> once
  std::vector<std::string> stops(stop_codons.size());
  for (int i = 0; i < stop_codons.size(); ++i) {
    stops[i] = Rcpp::as<std::string>(stop_codons[i]);
  }

  std::vector<int> out_starts, out_ends, out_indices;
  out_starts.reserve(seqs.size() * 4);
  out_ends.reserve(seqs.size() * 4);
  out_indices.reserve(seqs.size() * 4);

  for (int i = 0; i < seqs.size(); ++i) {
    if (CharacterVector::is_na(seqs[i])) continue;
    std::string seq = Rcpp::as<std::string>(seqs[i]);
    scan_seq(seq, start_codon, stops, min_body, i + 1,
             out_starts, out_ends, out_indices);
  }

  return List::create(
    Named("starts")  = IntegerVector(out_starts.begin(),  out_starts.end()),
    Named("ends")    = IntegerVector(out_ends.begin(),    out_ends.end()),
    Named("indices") = IntegerVector(out_indices.begin(), out_indices.end())
  );
}
