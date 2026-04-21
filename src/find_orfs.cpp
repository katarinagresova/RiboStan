#include <Rcpp.h>
using namespace Rcpp;

// Find all ORFs in a single nucleotide string.
//
// For each occurrence of ANY start codon in `start_codons`, walk in-frame
// (3 nt at a time) until a stop codon is encountered. Results are appended
// into the caller-supplied vectors `starts`, `ends`, `indices`, and
// `start_codon_idx` (index into `start_codons` for the matched start).
//
// Positions are 1-based and include both the start and stop codons.
// `min_body` is the minimum number of nucleotides between the end of the start
// codon and the beginning of the stop codon (0 = start + stop adjacent).
static void scan_seq(const std::string& seq,
                     const std::vector<std::string>& start_codons,
                     const std::vector<std::string>& stop_codons,
                     int min_body,
                     int seq_index,
                     std::vector<int>& starts,
                     std::vector<int>& ends,
                     std::vector<int>& indices,
                     std::vector<int>& start_codon_idx) {
  int n = (int)seq.size();

  for (std::size_t si = 0; si < start_codons.size(); ++si) {
    const std::string& start_codon = start_codons[si];
    int sc_len = (int)start_codon.size();   // normally 3
    std::size_t pos = 0;

    while ((pos = seq.find(start_codon, pos)) != std::string::npos) {
      int s = (int)pos;                       // 0-based start of start codon
      int body_start = s + sc_len;            // 0-based start of codon after start

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
            start_codon_idx.push_back((int)si + 1);  // 1-based R index
          }
          break;
        }
        p += 3;
      }
      ++pos;  // advance past current position to find all (possibly overlapping) starts
    }
  }
}


//' Find ORFs in a character vector of nucleotide sequences
//'
//' @param seqs Character vector of uppercase DNA sequences.
//' @param start_codons Character vector of start codons to search for
//'   (e.g. \code{c("ATG","CTG","GTG")}).
//' @param stop_codons Character vector; stop codons.
//' @param min_body Integer; minimum body length in nucleotides between the end
//'   of the start codon and the beginning of the stop codon.
//' @return A list with four integer vectors: \code{starts}, \code{ends}
//'   (1-based, inclusive), \code{indices} (1-based index into \code{seqs}),
//'   and \code{start_codon_idx} (1-based index into \code{start_codons}
//'   identifying which start codon was matched).
//' @keywords internal
// [[Rcpp::export]]
List find_orfs_cpp(CharacterVector seqs,
                   CharacterVector start_codons,
                   CharacterVector stop_codons,
                   int             min_body = 0) {
  std::vector<std::string> starts_vec(start_codons.size());
  for (int i = 0; i < start_codons.size(); ++i) {
    starts_vec[i] = Rcpp::as<std::string>(start_codons[i]);
  }
  std::vector<std::string> stops(stop_codons.size());
  for (int i = 0; i < stop_codons.size(); ++i) {
    stops[i] = Rcpp::as<std::string>(stop_codons[i]);
  }

  std::vector<int> out_starts, out_ends, out_indices, out_start_codon_idx;
  out_starts.reserve(seqs.size() * 4);
  out_ends.reserve(seqs.size() * 4);
  out_indices.reserve(seqs.size() * 4);
  out_start_codon_idx.reserve(seqs.size() * 4);

  for (int i = 0; i < seqs.size(); ++i) {
    if (CharacterVector::is_na(seqs[i])) continue;
    std::string seq = Rcpp::as<std::string>(seqs[i]);
    scan_seq(seq, starts_vec, stops, min_body, i + 1,
             out_starts, out_ends, out_indices, out_start_codon_idx);
  }

  return List::create(
    Named("starts")          = IntegerVector(out_starts.begin(),          out_starts.end()),
    Named("ends")            = IntegerVector(out_ends.begin(),            out_ends.end()),
    Named("indices")         = IntegerVector(out_indices.begin(),         out_indices.end()),
    Named("start_codon_idx") = IntegerVector(out_start_codon_idx.begin(), out_start_codon_idx.end())
  );
}
