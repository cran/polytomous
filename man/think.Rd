\name{think}
\alias{think}
\docType{data}
\title{
  Finnish \sQuote{think} verbs.
}
\description{
  3404 occurrences of four synonymous Finnish \sQuote{think} verbs (\sQuote{ajatella}: 1492; \sQuote{miettia}: 812; \sQuote{pohtia}: 713; \sQuote{harkita}: 387) in newspaper and Internet newsgroup discussion texts
}
\usage{data(think)}
\format{
  A data frame with 3404 observations on the following 27 variables.

  \describe{
    \item{\code{Lexeme}}{A factor specifying one of the four \sQuote{think} verb synonyms}
    \item{\code{Polarity}}{A factor specifying whether the \sQuote{think} verb has negative polarity or not (\code{=Other})}
    \item{\code{Voice}}{A factor specifying whether the \sQuote{think} verb is in the passive voice or not (\code{=Other})}
    \item{\code{Mood}}{A factor specifying whether the \sQuote{think} verb is in the indicative or conditional mood or not (=Other)}
    \item{\code{Person}}{A factor specifying whether the \sQuote{think} verb is in the first, second, third person or not (\code{=None})}
    \item{\code{Number}}{A factor specifying whether the \sQuote{think} verb is in the plural number or not (\code{=Other})}
    \item{\code{Covert}}{A factor specifying whether the agent/subject of the \sQuote{think} verb is explicitly expressed as a syntactic argument (\code{=Overt}), or only as a morphological feature of the \sQuote{think} verb (\code{=Covert})}
    \item{\code{ClauseEquivalent}}{A factor specifying whether the \sQuote{think} verb is used as a non-finite clause equivalent (\code{=ClauseEquivalent}) or as a finite verb for, (\code{=FiniteVerbChain})}
    \item{\code{Agent}}{A factor specifying the occurrence of Agent/Subject of the \sQuote{think} verb as either a Human Individual, Human Group, or as absent (\code{=None})}
    \item{\code{Patient}}{A factor specifying the occurrence of the Patient/Object argument among the semantic or structural subclasses as either an Human Individual/Group, Abstraction, Activity, Communication, an \sQuote{etta} (\sQuote{that}) clause (\code{=etta_CLAUSE}), DirectQuote, IndirectQuestion, Infinitive, Participle, or as absent (\code{=None})}
    \item{\code{Manner}}{A factor specifying the occurrrence of the Manner argument as any of its subclasses Generic, Negative (sufficiency), Positive (sufficiency), Frame, Agreement (Agree or Disagree), Joint, or as absent (\code{=None})}
    \item{\code{Time}}{A factor specifying the occurrence of Time argument (as a moment) as either of its subclasses definite (Definite), indefinite (Indefinite), or as absent (\code{=None})}
    \item{\code{Modality1}}{A factor specifying the main semantic subclasses of the entire Verb chain as either indicating Possibility, Necessity, or their absense (\code{=None})}
    \item{\code{Modality2}}{A factor specifying minor semantic subclasses of the entire Verb chain as indicating either a Temporal element (begin, end, continuation, etc.), External (cause), Volition, Accidental nature of the thinking process, or their absense (\code{=None})}
    \item{\code{Source}}{A factor specifying the occurrence of a Source argument or its absense (\code{=None})}
    \item{\code{Goal}}{A factor specifying the occurrence of a Goal argument or its absence (\code{=None})}
    \item{\code{Quantity}}{A factor specifying the occurrence of a Quantity argument, or its absence (\code{=None})}
    \item{\code{Location}}{A factor specifying the occurrence of a Location argument, or its absence (\code{=None})}
    \item{\code{Duration}}{A factor specifying the occurrence of a Duration argument, or its absence (\code{=None})}
    \item{\code{Frequency}}{A factor specifying the occurrence of a Frequency arument, or its absence (\code{=None})}
    \item{\code{MetaComment}}{A factor specifying the occurrence of a MetaComment, or its absence (\code{=None})}
    \item{\code{ReasonPurpose}}{A factor specifying the occurrence of a Reason or Purpose argument, or their absence (\code{=None})}
    \item{\code{Condition}}{A factor specifying the occurrence of a Condition argument, or its absence (\code{=None})}
    \item{\code{CoordinatedVerb}}{A factor specifying the occurrence of a Coordinated Verb (in relation to the \sQuote{think} verb), or its absence (\code{=None})}
    \item{\code{Register}}{A factor specifying whether the \sQuote{think} verb occurs in the newspaper subcorpus (\code{=hs95}) or the Internet newsgroup discussion corpus (\code{=sfnet})}
    \item{\code{Section}}{A factor specifying the subsection in which the \sQuote{think} verb occurs in either of the two subcorpora}
    \item{\code{Author}}{A factor specifying the author of the text in which the \sQuote{think} verb occurs, if that author is identifiable - authors in the Internet newgroup discussion subcorpus are anonymized; unidentifiable/unknown author designated as (\code{=None})}
}
}
\details{

The four most frequent synonyms meaning \sQuote{think, reflect,
ponder, consider}, i.e. \sQuote{ajatella, miettia, pohtia, harkita},
were extracted from two months of newspaper text from the 1990s
(Helsingin Sanomat 1995) and six months of Internet newsgroup
discussion from the early 2000s (SFNET 2002-2003), namely regarding
(personal) relationships (sfnet.keskustelu.ihmissuhteet) and politics
(sfnet.keskustelu.politiikka). The newspaper corpus consisted of
3,304,512 words of body text (i.e. excluding headers and captions as
well as punctuation tokens), and included 1,750 examples of the
studied \sQuote{think} verbs. The Internet corpus comprised 1,174,693
words of body text, yielding 1,654 instances of the selected
\sQuote{think} verbs. In terms of distinct identifiable authors, the
newspaper sub-corpus was the product of just over 500 journalists and
other contributors, while the Internet sub-corpus involved well over
1000 discussants. The \code{think} dataset contains a selection of 26
contextual features judged as most informative.

For extensive details of the data and its linguistic and statistical
analysis, see Arppe (2008). For the full selection of contextual
features, see the \code{amph} (2008) microcorpus.

}

\source{
amph (2008) A micro-corpus of 3404 occurrences of the four most common
Finnish THINK lexemes, `ajatella, miettia, pohtia, and harkita', in
Finnish newspaper and Internet newsgroup discussion texts, containing
extracts and linguistic analysis of the relevant context in the
original corpus data, scripts for processing this data, R functions
for its statistical analysis, as well as a comprehensive set of
ensuing results as R data tables. Compiled and analyzed by Antti
Arppe. Available on-line at URL:
http://www.csc.fi/english/research/software/amph/

Helsingin Sanomat (1995) ~22 million words of Finnish newspaper
articles published in Helsingin Sanomat during January--December
1995. Compiled by the Research Institute for the Languages of Finland
[KOTUS] and CSC -- IT Center for Science, Finland. Available on-line at
URL: http://www.csc.fi/kielipankki/

SFNET (2002-2003) ~100 million words of Finnish internet newsgroup
discussion posted during October 2002--April 2003. Compiled by Tuuli
Tuominen and Panu Kalliokoski, Computing Centre, University of
Helsinki, and Antti Arppe, Department of General Linguistics,
University of Helsinki, and CSC - IT Center for Science,
Finland. Available on-line at URL: http://www.csc.fi/kielipankki/
}
\references{

Arppe, A. (2008) Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Arppe, A. (2009) Linguistic choices vs. probabilities -- how much and
what can linguistic theory explain? In: Featherston, S. and S. Winkler
(eds.) The Fruits of Empirical Linguistics. Volume 1: Process. Berlin:
de Gruyter, pp. 1--24.

}
\examples{

## For examples see examples(polytomous)

}
\keyword{datasets}
