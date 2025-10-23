Diversity Analysis Important Notes
=================================
Data Type: Relative abundance data
Conversion Method: Relative abundance × 1,000,000 → round

ACE Calculation Issue:
- ACE index requires specific data conditions to calculate properly
- May fail if there are no singletons or specific abundance patterns
- Manual calculation method implemented as fallback

Indicator Reliability Classification:
[High Reliability]: Observed, Richness, Number, Shannon, Simpson, InvSimpson
    - Based on original relative abundance data
    - Results can be directly interpreted

[Interpret with Caution]: Chao1
    - Based on converted simulated count data
    - Affected by conversion factor

[If Available - Extra Caution]: ACE
    - Most sensitive to data conversion artifacts
    - Requires careful validation

Statistical Notes:
- Wilcoxon test used for two-group comparisons
- Kruskal-Wallis test used for multiple groups
- Tied values handled with normal approximation

Recommendation: Focus on high reliability indicators for main results
=================================
