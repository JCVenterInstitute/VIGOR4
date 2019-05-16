package org.jcvi.vigor.service;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.testing.category.Fast;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.testing.category.ReferenceDatabase;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.junit.ClassRule;
import org.junit.Rule;
import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.rules.SpringClassRule;
import org.springframework.test.context.junit4.rules.SpringMethodRule;

import java.io.File;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.notNullValue;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

@Category({Fast.class, ReferenceDatabase.class})
@RunWith(Parameterized.class)
@ContextConfiguration(classes = Application.class)
public class PeptideServiceTest {

    private final static Logger LOGGER = LogManager.getLogger(PeptideServiceTest.class);
    @ClassRule
    public static final SpringClassRule springClassRule = new SpringClassRule();
    @Rule
    public final SpringMethodRule springMethodRule = new SpringMethodRule();
    @Autowired
    private PeptideService peptideService;
    @Autowired
    private VigorInitializationService initializationService;
    private List<String[]> expected;
    private String id;
    private String sequence;
    private String mp_ref_db;

    public PeptideServiceTest ( String id, String sequence, String mp_ref_db, List<String[]> expected ) {

        this.id = id;
        this.sequence = sequence;
        this.mp_ref_db = mp_ref_db;
        this.expected = expected;
    }

    @Test
    public void testPeptides () throws VigorException {
        // must have atleast id, polyprotein, peptide db
        ProteinSequence protein = ProteinSequence.of(sequence);
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String refDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(refDBPath);
        assertThat("Reference database path must be set", refDBPath, notNullValue());
        File peptideDB = getPeptideDB(refDBPath, mp_ref_db);
        PeptideMatchingService.Scores scores = PeptideMatchingService.Scores.of(.25d, .4d, .5d);
        Model model = new Model();
        model.setTranslatedSeq(protein);
        model.setPartial3p(false);
        model.setPartial5p(false);
        Alignment alignment = new Alignment();
        ViralProtein viralProtein = new ViralProtein();
        viralProtein.setProteinID(id);
        alignment.setViralProtein(viralProtein);
        viralProtein.setSequence(protein);
        model.setAlignment(alignment);
        List<MaturePeptideMatch> matches = peptideService.findPeptides(model, peptideDB, scores);
        LOGGER.debug(() -> String.format("peptides:%s", matches.stream().map(String:: valueOf).collect(Collectors.joining("\n> ", "\n> ", ""))));
        Function<MaturePeptideMatch, ProteinSequence> matchToSeq = ( m ) -> {
            ProteinSequenceBuilder pb = m.getProtein().toBuilder().trim(m.getProteinRange());
            if (pb.get((int) pb.getLength() - 1) == AminoAcid.STOP) {
                pb.delete(Range.of(pb.getLength() - 1));
            }
            return pb.build();
        };
        if (matches.size() != expected.size()) {
            fail(String.format("for %s expected %s:\n%s\nactual %s:\n%s",
                    id,
                    expected.size(),
                    expected.stream().map(e -> e[ 2 ]).collect(Collectors.joining("\n")),
                    matches.size(),
                    matches.stream().map(m -> matchToSeq.apply(m).toString()).collect(Collectors.joining("\n"))
                    )
            );
        }
        String expectedProduct;
        String expectedID;
        ProteinSequence expectedSequence;
        MaturePeptideMatch match;
        ProteinSequence matchedSequence;
        Range prev = null;
        Range current;
        MaturePeptideMatch prevMatch = null;
        for (int i = 0; i < expected.size(); i++) {
            match = matches.get(i);
            expectedID = expected.get(i)[ 0 ];
            expectedProduct = expected.get(i)[ 1 ];
            expectedSequence = ProteinSequence.of(expected.get(i)[ 2 ]);
            matchedSequence = matchToSeq.apply(match);
            assertThat(matchedSequence, equalTo(expectedSequence));
            assertThat(match.getReference().getProduct(), equalTo(expectedProduct));
            current = match.getProteinRange();
            if (prev == null) {
                assertTrue(current.getBegin() == 0 || match.isFuzzyBegin());
            } else {
                assertTrue(String.format("ranges [%s:%s] fuzzyEnd=%s, fuzzyBegin=%s [%s:%s]  are not contiguous and not fuzzy",
                        prev.getBegin(), prev.getEnd(), prevMatch.isFuzzyEnd(),
                        match.isFuzzyBegin(), current.getBegin(), current.getEnd()), current.getBegin() == prev.getEnd() + 1 ||
                        ( match.isFuzzyBegin() && prevMatch.isFuzzyEnd() ));
            }
            prev = current;
            prevMatch = match;
        }
    }

    private File getPeptideDB ( String referenceDBPath, String mp_ref_db ) throws VigorException {

        File peptideDB = Paths.get(referenceDBPath, mp_ref_db).toFile();
        if (!peptideDB.exists()) {
            throw new VigorException(String.format("unable to find peptide DB for %s", mp_ref_db));
        }
        return peptideDB;
    }

    @Parameterized.Parameters(name = "PeptideServiceTest[#{index} {0}]")
    public static Collection<Object[]> getTestData () {

        List<Object[]> testData = new ArrayList<>();
        testData.add(
                new Object[] {
                        "HA_1_615",
                        String.join("",
                                "MKAILVVLLYTFATANADTLCIGYHANNST",
                                "DTVDTVLEKNVTVTHSVNLLEDKHNGKLCK",
                                "LRGVAPLHLGKCNIAGWILGNPECESLSTA",
                                "SSWSYIVETSSSDNGTCYPGDFINYEELRE",
                                "QLSSVSSFERFEIFPKTSSWPNHDSNKGVT",
                                "AACPQAGAKSFYKNLIWLVKKGNSYPKLSK",
                                "SYINDKGKEVLVLWGIHH"),
                        "flua_ha_mp",
                        Arrays.asList(
                                new String[] { "seg4matureA1", "signal peptide", "MKAILVVLLYTFATANA" },
                                new String[] { "seg4matureA2", "HA1",
                                        String.join("",
                                                "DTLCIGYHANNSTDTVDTVLEKNVTVTHSV",
                                                "NLLEDKHNGKLCKLRGVAPLHLGKCNIAGW",
                                                "ILGNPECESLSTASSWSYIVETSSSDNGTC",
                                                "YPGDFINYEELREQLSSVSSFERFEIFPKT",
                                                "SSWPNHDSNKGVTAACPQAGAKSFYKNLIW",
                                                "LVKKGNSYPKLSKSYINDKGKEVLVLWGIHH") })
                });
        testData.add(
                new Object[] {
                        "HA_1129_1753.1",
                        String.join("",
                                "QNAQGSGYAADLKSTQNAIDKITNKVNSVI",
                                "EKMNTQFTAVGKEFNHLEKRIENLNKKVDD",
                                "GFLDIWTYNAELLVLLENERTLDYHDSNVK",
                                "NLYEKVRSQLKNNAKEIGNGCFEFYHKCDN",
                                "TCMESVKNGTYDYPRYSEEAKLNREEIDGV",
                                "KLESTRIYQILAIYSTVASSLVLIVSLGAI",
                                "SFWMCSNGSLQCRICI*"),
                        "flua_ha_mp",
                        Arrays.<String[]> asList(new String[] { "seg4matureA3", "HA2",
                                String.join("",
                                        "QNAQGSGYAADLKSTQNAIDKITNKVNSVI",
                                        "EKMNTQFTAVGKEFNHLEKRIENLNKKVDD",
                                        "GFLDIWTYNAELLVLLENERTLDYHDSNVK",
                                        "NLYEKVRSQLKNNAKEIGNGCFEFYHKCDN",
                                        "TCMESVKNGTYDYPRYSEEAKLNREEIDGV",
                                        "KLESTRIYQILAIYSTVASSLVLIVSLGAI",
                                        "SFWMCSNGSLQCRICI") })
                });
        testData.add(
                new Object[] {
                        "gi_260907762.1",
                        String.join("",
                                "MNIQILVFALVAIIPTNADKICLGHHAVSN",
                                "GTKVNTLTERGVEVVNATETVERTNVPRIC",
                                "SKGKRTVDLGQCGLLGTITGPPQCDQFLEF",
                                "SADLIIERRGGSDVCYPGKFVNEEALRQIL",
                                "RESGGIDKETMGFTYSGIRTNGATSACRRS",
                                "GSSFYAEMKWLLSNTDNAAFPQMTKSYKNT",
                                "RKDPALIIWGIHHSGSTTEQTKLYGSGSKL",
                                "ITVGSSNYQQSFVPSPGARPQVNGQSGRID",
                                "FHWLMLNPNDTVTFSFNGAFIAPDRASFLK",
                                "GKSMGIQSGVQVDANCEGDCYHSGGTIISN",
                                "LPFQNINSRAVGKCPRYVKQESLMLATGMK",
                                "NVPELPKGRGLFGAIAGFIENGWEGLIDGW",
                                "YGFRHQNAQGEGTAADYKSTQSAIDQITGK",
                                "LNRLIEKTNQQFELIDNEFTEVEKQIGNVI",
                                "NWTRDSMTEVWSYNAELLVAMENQHTIDLA",
                                "DSEMNKLYERVKRQLRENAEEDGTGCFEIF",
                                "HKCDDDCMASIRNNTYDHSKYREEAMQNRI",
                                "QINPVKLSSGYKDVILWFSFGASCFILLAI",
                                "AMGLVFICVKNGNMRCTICI*"),
                        "flua_ha_mp",
                        Arrays.asList(new String[] { "seg4matureO1", "signal peptide",
                                        "MNIQILVFALVAIIPTNA" },
                                new String[] { "seg4mature7F2", "HA1",
                                        String.join("",
                                                "DKICLGHHAVSNGTKVNTLTERGVEVVNAT",
                                                "ETVERTNVPRICSKGKRTVDLGQCGLLGTI",
                                                "TGPPQCDQFLEFSADLIIERRGGSDVCYPG",
                                                "KFVNEEALRQILRESGGIDKETMGFTYSGI",
                                                "RTNGATSACRRSGSSFYAEMKWLLSNTDNA",
                                                "AFPQMTKSYKNTRKDPALIIWGIHHSGSTT",
                                                "EQTKLYGSGSKLITVGSSNYQQSFVPSPGA",
                                                "RPQVNGQSGRIDFHWLMLNPNDTVTFSFNG",
                                                "AFIAPDRASFLKGKSMGIQSGVQVDANCEG",
                                                "DCYHSGGTIISNLPFQNINSRAVGKCPRYV",
                                                "KQESLMLATGMKNVPELPKGR") },
                                new String[] { "seg4matureO3", "HA2",
                                        String.join("",
                                                "GLFGAIAGFIENGWEGLIDGWYGFRHQNAQ",
                                                "GEGTAADYKSTQSAIDQITGKLNRLIEKTN",
                                                "QQFELIDNEFTEVEKQIGNVINWTRDSMTE",
                                                "VWSYNAELLVAMENQHTIDLADSEMNKLYE",
                                                "RVKRQLRENAEEDGTGCFEIFHKCDDDCMA",
                                                "SIRNNTYDHSKYREEAMQNRIQINPVKLSS",
                                                "GYKDVILWFSFGASCFILLAIAMGLVFICV",
                                                "KNGNMRCTICI") })
                });
        testData.add(
                new Object[] {
                        "gi_155016323.1",
                        String.join("",
                                "MNTQILALIACMLIGAKGDKICLGHHAVAN",
                                "GTKVNTLTERGIEVVNATETVETANIKKIC",
                                "TQGKRPTDLGQCGLLGTLIGPPQCDQFLEF",
                                "DTDLIIERREGTDVCYPGKFTNEESLRQIL",
                                "RGSGGIDKESMGFTYSGIRTNGATSACRRS",
                                "GSSFYAEMKWLLSNSDNAAFPQMTKSYRNP",
                                "RNKPALIIWGVHHSGSATEQTKLYGSGNKL",
                                "ITVGSSKYQQSFTPSPGARPQVNGQSGRID",
                                "FHWLLLDPNDTVTFTFNGAFIAPDRASFFR",
                                "GESLGVQSDVPLDSGCEGDCFHSGGTIVSS",
                                "LPFQNINPRTVGKCPRYVKQTSLLLATGMR",
                                "NVPENPKTRGLFGAIAGFIENGWEGLIDGW",
                                "YGFRHQNAQGEGTAADYKSTQSAIDQITGK",
                                "LNRLIDKTNQQFELIDNEFSEIEQQIGNVI",
                                "NWTRDSMTEVWSYNAELLVAMENQHTIDLA",
                                "DSEMNKLYERVRKQLRENAEEDGTGCFEIF",
                                "HKCDDQCMESIRNNTYDHTQYRTESLQNRI",
                                "QIDPVRLSSGYKDIILWFSFGASCFLLLAI",
                                "AMGLVFICIKNGNMRCTICI*"),
                        "flua_ha_mp",
                        Arrays.asList(new String[] { "seg4mature7A1", "signal peptide",
                                        "MNTQILALIACMLIGAKG" },
                                new String[] { "seg4mature7C2", "HA1",
                                        String.join("",
                                                "DKICLGHHAVANGTKVNTLTERGIEVVNAT",
                                                "ETVETANIKKICTQGKRPTDLGQCGLLGTL",
                                                "IGPPQCDQFLEFDTDLIIERREGTDVCYPG",
                                                "KFTNEESLRQILRGSGGIDKESMGFTYSGI",
                                                "RTNGATSACRRSGSSFYAEMKWLLSNSDNA",
                                                "AFPQMTKSYRNPRNKPALIIWGVHHSGSAT",
                                                "EQTKLYGSGNKLITVGSSKYQQSFTPSPGA",
                                                "RPQVNGQSGRIDFHWLLLDPNDTVTFTFNG",
                                                "AFIAPDRASFFRGESLGVQSDVPLDSGCEG",
                                                "DCFHSGGTIVSSLPFQNINPRTVGKCPRYV",
                                                "KQTSLLLATGMRNVPENPKTR") },
                                new String[] { "seg4mature7A3", "HA2",
                                        String.join("",
                                                "GLFGAIAGFIENGWEGLIDGWYGFRHQNAQ",
                                                "GEGTAADYKSTQSAIDQITGKLNRLIDKTN",
                                                "QQFELIDNEFSEIEQQIGNVINWTRDSMTE",
                                                "VWSYNAELLVAMENQHTIDLADSEMNKLYE",
                                                "RVRKQLRENAEEDGTGCFEIFHKCDDQCME",
                                                "SIRNNTYDHTQYRTESLQNRIQIDPVRLSS",
                                                "GYKDIILWFSFGASCFLLLAIAMGLVFICI",
                                                "KNGNMRCTICI") }
                        )
                });
        return testData;
    }
}
